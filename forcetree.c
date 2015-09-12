#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <mpi.h>

#include "allvars.h"
#include "proto.h"


/*! \file forcetree.c
 *  \brief gravitational tree and code for Ewald correction
 *
 *  This file contains the computation of the gravitational force by means
 *  of a tree. The type of tree implemented is a geometrical oct-tree,
 *  starting from a cube encompassing all particles. This cube is
 *  automatically found in the domain decomposition, which also splits up
 *  the global "top-level" tree along node boundaries, moving the particles
 *  of different parts of the tree to separate processors. Tree nodes can
 *  be dynamically updated in drift/kick operations to avoid having to
 *  reconstruct the tree every timestep.
 */

/*! auxialiary variable used to set-up non-recursive walk */
static int last;



/*! length of lock-up table for short-range force kernel in TreePM algorithm */
#define NTAB 1000
/*! variables for short-range lookup table */
static float tabfac, shortrange_table[NTAB], shortrange_table_potential[NTAB];

/*! toggles after first tree-memory allocation, has only influence on log-files */
static int first_flag = 0;




#ifdef PERIODIC
/*! Macro that maps a distance to the nearest periodic neighbour */
#define NEAREST(x) (((x)>boxhalf)?((x)-boxsize):(((x)<-boxhalf)?((x)+boxsize):(x)))
/*! Size of 3D lock-up table for Lattice sum correction force */
#define EN  64
/*! 3D lock-up table for Lattice sum correction to force and potential. Only one
 *  octant is stored, the rest constructed by using the symmetry
 */
static FLOAT fcorrx[N_GRAVS][N_GRAVS][EN + 1][EN + 1][EN + 1];
static FLOAT fcorry[N_GRAVS][N_GRAVS][EN + 1][EN + 1][EN + 1];
static FLOAT fcorrz[N_GRAVS][N_GRAVS][EN + 1][EN + 1][EN + 1];
static FLOAT potcorr[N_GRAVS][N_GRAVS][EN + 1][EN + 1][EN + 1];
static double fac_intp;
#endif



/*! This function is a driver routine for constructing the gravitational
 *  oct-tree, which is done by calling a small number of other functions.
 */
int force_treebuild(int npart)
{
  Numnodestree = force_treebuild_single(npart);

  force_update_pseudoparticles();

  force_flag_localnodes();

  TimeOfLastTreeConstruction = All.Time;

  return Numnodestree;
}



/*! Constructs the gravitational oct-tree.  
 *
 *  The index convention for accessing tree nodes is the following: the
 *  indices 0...NumPart-1 reference single particles, the indices
 *  All.MaxPart.... All.MaxPart+nodes-1 reference tree nodes. `Nodes_base'
 *  points to the first tree node, while `nodes' is shifted such that
 *  nodes[All.MaxPart] gives the first tree node. Finally, node indices
 *  with values 'All.MaxPart + MaxNodes' and larger indicate "pseudo
 *  particles", i.e. multipole moments of top-level nodes that lie on
 *  different CPUs. If such a node needs to be opened, the corresponding
 *  particle must be exported to that CPU. The 'Extnodes' structure
 *  parallels that of 'Nodes'. Its information is only needed for the SPH
 *  part of the computation. (The data is split onto these two structures
 *  as a tuning measure.  If it is merged into 'Nodes' a somewhat bigger
 *  size of the nodes also for gravity would result, which would reduce
 *  cache utilization slightly.
 */
int force_treebuild_single(int npart)
{
  int i, j, subnode = 0, parent, numnodes;
  int nfree, th, nn, no;
  struct NODE *nfreep;
  double lenhalf, epsilon;
  peanokey key;


  /* create an empty root node  */
  nfree = All.MaxPart;		/* index of first free node */
  nfreep = &Nodes[nfree];	/* select first node */

  nfreep->len = DomainLen;
  for(j = 0; j < 3; j++)
    nfreep->center[j] = DomainCenter[j];
  for(j = 0; j < 8; j++)
    nfreep->u.suns[j] = -1;


  numnodes = 1;
  nfreep++;
  nfree++;

  /* create a set of empty nodes corresponding to the top-level domain
   * grid. We need to generate these nodes first to make sure that we have a
   * complete top-level tree which allows the easy insertion of the
   * pseudo-particles at the right place 
   */

  force_create_empty_nodes(All.MaxPart, 0, 1, 0, 0, 0, &numnodes, &nfree);


  /* if a high-resolution region in a global tree is used, we need to generate
   * an additional set empty nodes to make sure that we have a complete
   * top-level tree for the high-resolution inset
   */

  nfreep = &Nodes[nfree];
  parent = -1;			/* note: will not be used below before it is changed */


  /* now we insert all particles */
  for(i = 0; i < npart; i++)
    {

      /* the softening is only used to check whether particles are so close
       * that the tree needs not to be refined further
       */
      epsilon = All.ForceSoftening[P[i].Type];

      key = peano_hilbert_key((P[i].Pos[0] - DomainCorner[0]) * DomainFac,
			      (P[i].Pos[1] - DomainCorner[1]) * DomainFac,
			      (P[i].Pos[2] - DomainCorner[2]) * DomainFac, BITS_PER_DIMENSION);

      no = 0;
      while(TopNodes[no].Daughter >= 0)
	no = TopNodes[no].Daughter + (key - TopNodes[no].StartKey) / (TopNodes[no].Size / 8);

      no = TopNodes[no].Leaf;
      th = DomainNodeIndex[no];

      while(1)
	{
	  if(th >= All.MaxPart)	/* we are dealing with an internal node */
	    {
	      subnode = 0;
	      if(P[i].Pos[0] > Nodes[th].center[0])
		subnode += 1;
	      if(P[i].Pos[1] > Nodes[th].center[1])
		subnode += 2;
	      if(P[i].Pos[2] > Nodes[th].center[2])
		subnode += 4;

	      nn = Nodes[th].u.suns[subnode];

	      if(nn >= 0)	/* ok, something is in the daughter slot already, need to continue */
		{
		  parent = th;
		  th = nn;
		}
	      else
		{
		  /* here we have found an empty slot where we can attach
		   * the new particle as a leaf.
		   */
		  Nodes[th].u.suns[subnode] = i;
		  break;	/* done for this particle */
		}
	    }
	  else
	    {
	      /* We try to insert into a leaf with a single particle.  Need
	       * to generate a new internal node at this point.
	       */
	      Nodes[parent].u.suns[subnode] = nfree;

	      nfreep->len = 0.5 * Nodes[parent].len;
	      lenhalf = 0.25 * Nodes[parent].len;

	      if(subnode & 1)
		nfreep->center[0] = Nodes[parent].center[0] + lenhalf;
	      else
		nfreep->center[0] = Nodes[parent].center[0] - lenhalf;

	      if(subnode & 2)
		nfreep->center[1] = Nodes[parent].center[1] + lenhalf;
	      else
		nfreep->center[1] = Nodes[parent].center[1] - lenhalf;

	      if(subnode & 4)
		nfreep->center[2] = Nodes[parent].center[2] + lenhalf;
	      else
		nfreep->center[2] = Nodes[parent].center[2] - lenhalf;

	      nfreep->u.suns[0] = -1;
	      nfreep->u.suns[1] = -1;
	      nfreep->u.suns[2] = -1;
	      nfreep->u.suns[3] = -1;
	      nfreep->u.suns[4] = -1;
	      nfreep->u.suns[5] = -1;
	      nfreep->u.suns[6] = -1;
	      nfreep->u.suns[7] = -1;


	      subnode = 0;
	      if(P[th].Pos[0] > nfreep->center[0])
		subnode += 1;
	      if(P[th].Pos[1] > nfreep->center[1])
		subnode += 2;
	      if(P[th].Pos[2] > nfreep->center[2])
		subnode += 4;
#ifndef NOTREERND
	      if(nfreep->len < 1.0e-3 * epsilon)
		{
		  /* seems like we're dealing with particles at identical (or extremely close)
		   * locations. Randomize subnode index to allow tree construction. Note: Multipole moments
		   * of tree are still correct, but this will only happen well below gravitational softening
		   * length-scale anyway.
		   */
		  subnode = (int) (8.0 * get_random_number((0xffff & P[i].ID) + P[i].GravCost));
		  P[i].GravCost += 1;
		  if(subnode >= 8)
		    subnode = 7;
		}
#endif
	      nfreep->u.suns[subnode] = th;

	      th = nfree;	/* resume trying to insert the new particle at
				 * the newly created internal node
				 */

	      numnodes++;
	      nfree++;
	      nfreep++;

	      if((numnodes) >= MaxNodes)
		{
		  printf("task %d: maximum number %d of tree-nodes reached.\n", ThisTask, MaxNodes);
		  printf("for particle %d\n", i);
		  dump_particles();
		  endrun(1);
		}
	    }
	}
    }


  /* insert the pseudo particles that represent the mass distribution of other domains */
  force_insert_pseudo_particles();


  /* now compute the multipole moments recursively */
  last = -1;

  force_update_node_recursive(All.MaxPart, -1, -1);

  if(last >= All.MaxPart)
    {
      if(last >= All.MaxPart + MaxNodes)	/* a pseudo-particle */
	Nextnode[last - MaxNodes] = -1;
      else
	Nodes[last].u.d.nextnode = -1;
    }
  else
    Nextnode[last] = -1;

  return numnodes;
}



/*! This function recursively creates a set of empty tree nodes which
 *  corresponds to the top-level tree for the domain grid. This is done to
 *  ensure that this top-level tree is always "complete" so that we can
 *  easily associate the pseudo-particles of other CPUs with tree-nodes at
 *  a given level in the tree, even when the particle population is so
 *  sparse that some of these nodes are actually empty.
*/
void force_create_empty_nodes(int no, int topnode, int bits, int x, int y, int z, int *nodecount,
			      int *nextfree)
{
  int i, j, k, n, sub, count;

  if(TopNodes[topnode].Daughter >= 0)
    {
      for(i = 0; i < 2; i++)
	for(j = 0; j < 2; j++)
	  for(k = 0; k < 2; k++)
	    {
	      sub = 7 & peano_hilbert_key((x << 1) + i, (y << 1) + j, (z << 1) + k, bits);

	      count = i + 2 * j + 4 * k;

	      Nodes[no].u.suns[count] = *nextfree;


	      Nodes[*nextfree].len = 0.5 * Nodes[no].len;
	      Nodes[*nextfree].center[0] = Nodes[no].center[0] + (2 * i - 1) * 0.25 * Nodes[no].len;
	      Nodes[*nextfree].center[1] = Nodes[no].center[1] + (2 * j - 1) * 0.25 * Nodes[no].len;
	      Nodes[*nextfree].center[2] = Nodes[no].center[2] + (2 * k - 1) * 0.25 * Nodes[no].len;

	      for(n = 0; n < 8; n++)
		Nodes[*nextfree].u.suns[n] = -1;

	      if(TopNodes[TopNodes[topnode].Daughter + sub].Daughter == -1)
		DomainNodeIndex[TopNodes[TopNodes[topnode].Daughter + sub].Leaf] = *nextfree;

	      *nextfree = *nextfree + 1;
	      *nodecount = *nodecount + 1;

	      if((*nodecount) >= MaxNodes)
		{
		  printf("task %d: maximum number %d of tree-nodes reached.\n", ThisTask, MaxNodes);
		  printf("in create empty nodes\n");
		  dump_particles();
		  endrun(11);
		}

	      force_create_empty_nodes(*nextfree - 1, TopNodes[topnode].Daughter + sub,
				       bits + 1, 2 * x + i, 2 * y + j, 2 * z + k, nodecount, nextfree);
	    }
    }
}


/*! this function inserts pseudo-particles which will represent the mass
 *  distribution of the other CPUs. Initially, the mass of the
 *  pseudo-particles is set to zero, and their coordinate is set to the
 *  center of the domain-cell they correspond to. These quantities will be
 *  updated later on.
 */
void force_insert_pseudo_particles(void)
{
  int i, index, subnode, nn, th;
  int j;

  for(i = 0; i < NTopleaves; i++)
    {
      index = DomainNodeIndex[i];
      
      for(j = 0; j < N_GRAVS; ++j) {
	DomainMoment[i].s[0][j] = Nodes[index].center[0];
	DomainMoment[i].s[1][j] = Nodes[index].center[1];
	DomainMoment[i].s[2][j] = Nodes[index].center[2];
	DomainMoment[i].mass[j] = 0;
	
#ifdef NGRAVS_ACCUMULATOR
	DomainMoment[i].Nparticles[j] = 0;
#endif
      }
    }

  for(i = 0; i < NTopleaves; i++)
    {
      if(i < DomainMyStart || i > DomainMyLast)
	{
	  th = All.MaxPart;	/* select index of first node in tree */

	  while(1)
	    {
	      if(th >= All.MaxPart)	/* we are dealing with an internal node */
		{
		  if(th >= All.MaxPart + MaxNodes)
		    endrun(888);	/* this can't be */

		  subnode = 0;
		  
		  // KC 8/11/14 This looks like a bisection check, so the generalization should be that
		  // if any of the mass center coordinates exceeds the center, than we subnode up
		  for(j = 0; j < N_GRAVS; ++j) {
		    
		    if(DomainMoment[i].s[0][j] > Nodes[th].center[0]) {
		      subnode += 1;
		      break;
		    }
		  }

		  for(j = 0; j < N_GRAVS; ++j) {
		  
		    if(DomainMoment[i].s[1][j] > Nodes[th].center[1]) {
		      subnode += 2;
		      break;
		    }
		  }

		  for(j = 0; j < N_GRAVS; ++j) {
		  
		    if(DomainMoment[i].s[2][j] > Nodes[th].center[2]) {
		      subnode += 4;
		      break;
		    }
		  }

		  //fprintf(stderr, "subnode structure: [%d, %d, %d]\n", i, th, subnode);
		  nn = Nodes[th].u.suns[subnode];

		  if(nn >= 0)	/* ok, something is in the daughter slot already, need to continue */
		    {
		      th = nn;
		    }
		  else
		    {
		      /* here we have found an empty slot where we can 
		       * attach the pseudo particle as a leaf 
		       */
		      Nodes[th].u.suns[subnode] = All.MaxPart + MaxNodes + i;

		      break;	/* done for this pseudo particle */
		    }
		}
	      else
		{
		  endrun(889);	/* this can't be */
		}
	    }
	}
    }
}

// KC 8/10/14  Now we update this function to correctly determine the (monopole) moments
// of all the different types of interactions we can have

/*! this routine determines the multipole moments for a given internal node
 *  and all its subnodes using a recursive computation.  The result is
 *  stored in the Nodes[] structure in the sequence of this tree-walk.
 *
 *  Note that the bitflags-variable for each node is used to store in the
 *  lowest bits some special information: Bit 0 flags whether the node
 *  belongs to the top-level tree corresponding to the domain
 *  decomposition, while Bit 1 signals whether the top-level node is
 *  dependent on local mass.
 * 
 *  If UNEQUALSOFTENINGS is set, bits 2-4 give the particle type with
 *  the maximum softening among the particles in the node, and bit 5
 *  flags whether the node contains any particles with lower softening
 *  than that.  
 */
void force_update_node_recursive(int no, int sib, int father)
{
  int j, jj, p, pp, nextsib, suns[8];
  FLOAT hmax;

#ifdef UNEQUALSOFTENINGS
#ifndef ADAPTIVE_GRAVSOFT_FORGAS
  int maxsofttype, diffsoftflag;
#else
  FLOAT maxsoft;
#endif
#endif
  struct particle_data *pa;
  int i;

  // KC 8/11/14 Updated to reflect the different types
  // ??? Why is this a double explicitly and not a FLOAT?
  // Perhaps it is because we are adding many things, and so we don't 
  // want to lose precision during the addition...
  double s[3][N_GRAVS], vs[3][N_GRAVS], mass[N_GRAVS];
#ifdef NGRAVS_ACCUMULATOR
  long Nparticles[N_GRAVS];
#endif

  if(no >= All.MaxPart && no < All.MaxPart + MaxNodes)	/* internal node */
    {
      for(j = 0; j < 8; j++)
	suns[j] = Nodes[no].u.suns[j];	/* this "backup" is necessary because the nextnode entry will
					   overwrite one element (union!) */
      if(last >= 0)
	{
	  if(last >= All.MaxPart)
	    {
	      if(last >= All.MaxPart + MaxNodes)	/* a pseudo-particle */
		Nextnode[last - MaxNodes] = no;
	      else
		Nodes[last].u.d.nextnode = no;
	    }
	  else
	    Nextnode[last] = no;
	}

      last = no;

      
      // KC 3/11/14
      // Replace memsets, memcpys with explicit looping.  3-9% performance gains!
      // Unroll inner loop: 1-2% performance gains!  Holy shit.
      for(i = 0; i < N_GRAVS; ++i) {

	s[0][i] = 0;
	s[1][i] = 0;
	s[2][i] = 0;

	vs[0][i] = 0;
	vs[1][i] = 0;
	vs[2][i] = 0;
	
	mass[i] = 0;
#ifdef NGRAVS_ACCUMULATOR
	Nparticles[i] = 0;
#endif
      }

     
      hmax = 0;
#ifdef UNEQUALSOFTENINGS
#ifndef ADAPTIVE_GRAVSOFT_FORGAS
      maxsofttype = 7;
      diffsoftflag = 0;
#else
      maxsoft = 0;
#endif
#endif

      for(j = 0; j < 8; j++)
	{
	  if((p = suns[j]) >= 0)
	    {
	      /* check if we have a sibling on the same level */
	      for(jj = j + 1; jj < 8; jj++)
		if((pp = suns[jj]) >= 0)
		  break;

	      if(jj < 8)	/* yes, we do */
		nextsib = pp;
	      else
		nextsib = sib;

	      force_update_node_recursive(p, nextsib, no);	      

	      if(p >= All.MaxPart)	/* an internal node or pseudo particle */
		{
		  if(p >= All.MaxPart + MaxNodes)	/* a pseudo particle */
		    {
		      /* nothing to be done here because the mass of the
		       * pseudo-particle is still zero. This will be changed
		       * later.
		       */
		    }
		  else
		    {
		      // KC 8/11/14 Loop over interaction types
		      // Note we cannot just memcpy due to type differences
		      for(i = 0; i < N_GRAVS; ++i) {
			
#ifdef NGRAVS_ACCUMULATOR	
			Nparticles[i] += Nodes[p].u.d.Nparticles[i];
#endif
			mass[i] += Nodes[p].u.d.mass[i];
			s[0][i] += Nodes[p].u.d.mass[i] * Nodes[p].u.d.s[0][i];
			s[1][i] += Nodes[p].u.d.mass[i] * Nodes[p].u.d.s[1][i];
			s[2][i] += Nodes[p].u.d.mass[i] * Nodes[p].u.d.s[2][i];
			vs[0][i] += Nodes[p].u.d.mass[i] * Extnodes[p].vs[0][i];
			vs[1][i] += Nodes[p].u.d.mass[i] * Extnodes[p].vs[1][i];
			vs[2][i] += Nodes[p].u.d.mass[i] * Extnodes[p].vs[2][i];
		      }

		      
		      if(Extnodes[p].hmax > hmax)
			hmax = Extnodes[p].hmax;

#ifdef UNEQUALSOFTENINGS
#ifndef ADAPTIVE_GRAVSOFT_FORGAS
		      diffsoftflag |= (Nodes[p].u.d.bitflags >> 5) & 1;

		      if(maxsofttype == 7)
			{
			  maxsofttype = (Nodes[p].u.d.bitflags >> 2) & 7;
			}
		      else
			{
			  if(((Nodes[p].u.d.bitflags >> 2) & 7) != 7)
			    {
			      if(All.ForceSoftening[((Nodes[p].u.d.bitflags >> 2) & 7)] >
				 All.ForceSoftening[maxsofttype])
				{
				  maxsofttype = ((Nodes[p].u.d.bitflags >> 2) & 7);
				  diffsoftflag = 1;
				}
			      else
				{
				  if(All.ForceSoftening[((Nodes[p].u.d.bitflags >> 2) & 7)] <
				     All.ForceSoftening[maxsofttype])
				    diffsoftflag = 1;
				}
			    }
			}
#else
		      if(Nodes[p].maxsoft > maxsoft)
			maxsoft = Nodes[p].maxsoft;
#endif
#endif
		    }
		}
	      else		/* a particle */
		{
		  pa = &P[p];
		  
		  
		  i = TypeToGrav[pa->Type];

		  mass[i] += pa->Mass;
		  s[0][i] += pa->Mass * pa->Pos[0];
		  s[1][i] += pa->Mass * pa->Pos[1];
		  s[2][i] += pa->Mass * pa->Pos[2];
		  vs[0][i] += pa->Mass * pa->Vel[0];
		  vs[1][i] += pa->Mass * pa->Vel[1];
		  vs[2][i] += pa->Mass * pa->Vel[2];

#ifdef NGRAVS_ACCUMULATOR		  
		  Nparticles[i]++;
#endif

#ifdef UNEQUALSOFTENINGS
#ifndef ADAPTIVE_GRAVSOFT_FORGAS
		  if(maxsofttype == 7)
		    {
		      maxsofttype = pa->Type;
		    }
		  else
		    {
		      if(All.ForceSoftening[pa->Type] > All.ForceSoftening[maxsofttype])
			{
			  maxsofttype = pa->Type;
			  diffsoftflag = 1;
			}
		      else
			{
			  if(All.ForceSoftening[pa->Type] < All.ForceSoftening[maxsofttype])
			    diffsoftflag = 1;
			}
		    }
#else
		  if(pa->Type == 0)
		    {
		      if(SphP[p].Hsml > maxsoft)
			maxsoft = SphP[p].Hsml;
		    }
		  else
		    {
		      if(All.ForceSoftening[pa->Type] > maxsoft)
			maxsoft = All.ForceSoftening[pa->Type];
		    }
#endif
#endif
		  if(pa->Type == 0)
		    if(SphP[p].Hsml > hmax)
		      hmax = SphP[p].Hsml;
		}
	    }
	}

      // KC 8/11/14 Looks like we should be looping over this stuff
      // KC 2/7/15 PPP (perform a single division, then the others are multuplications)
      // but Volker does this too, so this can't be so bad?
      for(i = 0; i < N_GRAVS; ++i) {

	if(mass[i] > 0)
	  {
	    s[0][i] /= mass[i];
	    s[1][i] /= mass[i];
	    s[2][i] /= mass[i];
	    vs[0][i] /= mass[i];
	    vs[1][i] /= mass[i];
	    vs[2][i] /= mass[i];
	  }
	else
	  {
	    s[0][i] = Nodes[no].center[0];
	    s[1][i] = Nodes[no].center[1];
	    s[2][i] = Nodes[no].center[2];
	  }
	
	Nodes[no].u.d.s[0][i] = s[0][i];
	Nodes[no].u.d.s[1][i] = s[1][i];
	Nodes[no].u.d.s[2][i] = s[2][i];
	Nodes[no].u.d.mass[i] = mass[i];

#ifdef NGRAVS_ACCUMULATOR
	Nodes[no].u.d.Nparticles[i] = Nparticles[i];
#endif

	// KC 8/11/14 We move the vs up here as we need to loop to set all these things      
	Extnodes[no].vs[0][i] = vs[0][i];
	Extnodes[no].vs[1][i] = vs[1][i];
	Extnodes[no].vs[2][i] = vs[2][i];

#ifdef NGRAVS_ACCUMULATOR_DEBUG
	printf("%d: {%g, %g, %g} --> %g = %d\n", i, Nodes[no].center[0], Nodes[no].center[1], Nodes[no].center[2], Nodes[no].len, Nparticles[i]); 
#endif
      }

#ifdef UNEQUALSOFTENINGS
#ifndef ADAPTIVE_GRAVSOFT_FORGAS
      Nodes[no].u.d.bitflags = 4 * maxsofttype + 32 * diffsoftflag;
#else
      Nodes[no].u.d.bitflags = 0;
      Nodes[no].maxsoft = maxsoft;
#endif
#else
      Nodes[no].u.d.bitflags = 0;
#endif

      // KC 8/11/14 This is where the vs used to be
      Extnodes[no].hmax = hmax;

      Nodes[no].u.d.sibling = sib;
      Nodes[no].u.d.father = father;

    }
  else				/* single particle or pseudo particle */
    {
      if(last >= 0)
	{
	  if(last >= All.MaxPart)
	    {
	      if(last >= All.MaxPart + MaxNodes)	/* a pseudo-particle */
		Nextnode[last - MaxNodes] = no;
	      else
		Nodes[last].u.d.nextnode = no;
	    }
	  else
	    Nextnode[last] = no;
	}

      last = no;

      if(no < All.MaxPart)	/* only set it for single particles */
	Father[no] = father;
    }

}



/*! This function updates the multipole moments of the pseudo-particles
 *  that represent the mass distribution on different CPUs. For that
 *  purpose, it first exchanges the necessary data, and then updates the
 *  top-level tree accordingly. The detailed implementation of these two
 *  tasks is done in separate functions.
 */
void force_update_pseudoparticles(void)
{
  force_exchange_pseudodata();

  force_treeupdate_pseudos();
}



/*! This function communicates the values of the multipole moments of the
 *  top-level tree-nodes of the domain grid.  This data can then be used to
 *  update the pseudo-particles on each CPU accordingly.
 */
void force_exchange_pseudodata(void)
{
  int i, no;
  MPI_Status status;
  int level, sendTask, recvTask;
  int n, k;

  for(i = DomainMyStart; i <= DomainMyLast; i++)
    {
      no = DomainNodeIndex[i];

      /* read out the multipole moments from the local base cells */
      // KC 11/2/14
      // Unroll like a beast
      for(n = 0; n < N_GRAVS; ++n) {

	DomainMoment[i].s[0][n] = Nodes[no].u.d.s[0][n];
	DomainMoment[i].s[1][n] = Nodes[no].u.d.s[1][n];
	DomainMoment[i].s[2][n] = Nodes[no].u.d.s[2][n];
	DomainMoment[i].mass[n] = Nodes[no].u.d.mass[n];
#ifdef NGRAVS_ACCUMULATOR
	DomainMoment[i].Nparticles[n] = Nodes[no].u.d.Nparticles[n];
#endif
	DomainMoment[i].vs[0][n] = Extnodes[no].vs[0][n];
	DomainMoment[i].vs[1][n] = Extnodes[no].vs[1][n];
	DomainMoment[i].vs[2][n] = Extnodes[no].vs[2][n];
      }

#ifdef UNEQUALSOFTENINGS
#ifndef ADAPTIVE_GRAVSOFT_FORGAS
      DomainMoment[i].bitflags = Nodes[no].u.d.bitflags;
#else
      DomainMoment[i].maxsoft = Nodes[no].maxsoft;
#endif
#endif
    }

  /* share the pseudo-particle data accross CPUs */

  for(level = 1; level < (1 << PTask); level++)
    {
      sendTask = ThisTask;
      recvTask = ThisTask ^ level;

      if(recvTask < NTask)
	MPI_Sendrecv(&DomainMoment[DomainStartList[sendTask]],
		     (DomainEndList[sendTask] - DomainStartList[sendTask] + 1) * sizeof(struct DomainNODE),
		     MPI_BYTE, recvTask, TAG_DMOM,
		     &DomainMoment[DomainStartList[recvTask]],
		     (DomainEndList[recvTask] - DomainStartList[recvTask] + 1) * sizeof(struct DomainNODE),
		     MPI_BYTE, recvTask, TAG_DMOM, MPI_COMM_WORLD, &status);
    }

}

/*! This function updates the top-level tree after the multipole moments of
 *  the pseudo-particles have been updated.
 */
void force_treeupdate_pseudos(void)
{
  int i, k, no;
  // KC 8/11/14 Extension of the variables
  FLOAT sold[3][N_GRAVS], vsold[3][N_GRAVS], snew[3][N_GRAVS], vsnew[3][N_GRAVS], massold[N_GRAVS], massnew[N_GRAVS], mm[N_GRAVS];

#ifdef UNEQUALSOFTENINGS
#ifndef ADAPTIVE_GRAVSOFT_FORGAS
  int maxsofttype, diffsoftflag;
#else
  FLOAT maxsoft;
#endif
#endif
  int j;
#ifdef NGRAVS_ACCUMULATOR
  int partnew[N_GRAVS], partold[N_GRAVS];
#endif

  for(i = 0; i < NTopleaves; i++)
    if(i < DomainMyStart || i > DomainMyLast)
      {
	no = DomainNodeIndex[i];
	
	// KC 11/4/14
	// No need to unroll here because we keep hitting different variables each time
	// so little cache gain?
	//
	// He does keep it so that DomainMoment occurs near in the code
	for(j = 0; j < N_GRAVS; j++) {

	  for(k = 0; k < 3; ++k) {
	    vsold[k][j] = Extnodes[no].vs[k][j];
	    sold[k][j] = Nodes[no].u.d.s[k][j];
	  }
	  massold[j] = Nodes[no].u.d.mass[j];
#ifdef NGRAVS_ACCUMULATOR
	  partold[j] = Nodes[no].u.d.Nparticles[j];
#endif
	  for(k = 0; k < 3; ++k) {
	    
	    snew[k][j] = DomainMoment[i].s[k][j];
	    vsnew[k][j] = DomainMoment[i].vs[k][j];
	  }
	  massnew[j] = DomainMoment[i].mass[j];
#ifdef NGRAVS_ACCUMULATOR
	  partnew[j] = DomainMoment[i].Nparticles[j];
#endif
        }

#ifdef UNEQUALSOFTENINGS
#ifndef ADAPTIVE_GRAVSOFT_FORGAS
	maxsofttype = (DomainMoment[i].bitflags >> 2) & 7;
	diffsoftflag = (DomainMoment[i].bitflags >> 5) & 1;
#else
	maxsoft = DomainMoment[i].maxsoft;
#endif
#endif
	do
	  {
	    // KC 11/4/14
	    // This should correctly compute mass centers, as long as the values accumulated
	    // in mass[i] represent asymptotic masses of the sourcing entities
	    for(j = 0; j < N_GRAVS; ++j) {

#ifdef NGRAVS_ACCUMULATOR
	      Nodes[no].u.d.Nparticles[j] += partnew[j] - partold[j];
#endif
	      mm[j] = Nodes[no].u.d.mass[j] + massnew[j] - massold[j];

	      for(k = 0; k < 3; k++)
		{
		  if(mm[j] > 0)
		    {
		      Nodes[no].u.d.s[k][j] =
			(Nodes[no].u.d.mass[j] * Nodes[no].u.d.s[k][j] + 
			 massnew[j] * snew[k][j] - 
			 massold[j] * sold[k][j]) / mm[j];
		      Extnodes[no].vs[k][j] =
			(Nodes[no].u.d.mass[j] * Extnodes[no].vs[k][j] + 
			 massnew[j] * vsnew[k][j] -
			 massold[j] * vsold[k][j]) / mm[j];
		    }
		}
	      Nodes[no].u.d.mass[j] = mm[j];
	    }

#ifdef UNEQUALSOFTENINGS
#ifndef ADAPTIVE_GRAVSOFT_FORGAS
	    diffsoftflag |= (Nodes[no].u.d.bitflags >> 5) & 1;

	    if(maxsofttype == 7)
	      maxsofttype = (Nodes[no].u.d.bitflags >> 2) & 7;
	    else
	      {
		if(((Nodes[no].u.d.bitflags >> 2) & 7) != 7)
		  {
		    if(All.ForceSoftening[((Nodes[no].u.d.bitflags >> 2) & 7)] >
		       All.ForceSoftening[maxsofttype])
		      {
			maxsofttype = ((Nodes[no].u.d.bitflags >> 2) & 7);
			diffsoftflag = 1;
		      }
		    else
		      {
			if(All.ForceSoftening[((Nodes[no].u.d.bitflags >> 2) & 7)] <
			   All.ForceSoftening[maxsofttype])
			  diffsoftflag = 1;
		      }
		  }
	      }

	    Nodes[no].u.d.bitflags = (Nodes[no].u.d.bitflags & 3) + 4 * maxsofttype + 32 * diffsoftflag;
#else
	    if(Nodes[no].maxsoft < maxsoft)
	      Nodes[no].maxsoft = maxsoft;
	    maxsoft = Nodes[no].maxsoft;
#endif
#endif
	    no = Nodes[no].u.d.father;

	  }
	while(no >= 0);
      }
}



/*! This function flags nodes in the top-level tree that are dependent on
 *  local particle data.
 */
void force_flag_localnodes(void)
{
  int no, i;

  /* mark all top-level nodes */

  for(i = 0; i < NTopleaves; i++)
    {
      no = DomainNodeIndex[i];

      while(no >= 0)
	{
	  if((Nodes[no].u.d.bitflags & 1))
	    break;

	  Nodes[no].u.d.bitflags |= 1;

	  no = Nodes[no].u.d.father;
	}
    }

  /* mark top-level nodes that contain local particles */

  for(i = DomainMyStart; i <= DomainMyLast; i++)
    {
      /*
         if(DomainMoment[i].mass > 0)
       */
      {
	no = DomainNodeIndex[i];

	while(no >= 0)
	  {
	    if((Nodes[no].u.d.bitflags & 2))
	      break;

	    Nodes[no].u.d.bitflags |= 2;

	    no = Nodes[no].u.d.father;
	  }
      }
    }
}



/*! This function updates the side-length of tree nodes in case the tree is
 *  not reconstructed, but only drifted.  The grouping of particles to tree
 *  nodes is not changed in this case, but some tree nodes may need to be
 *  enlarged because particles moved out of their original bounds.
 */
void force_update_len(void)
{
  int i, no;
  MPI_Status status;
  int level, sendTask, recvTask;

  force_update_node_len_local();

  /* first update the side-lengths of all local nodes */
  for(i = DomainMyStart; i <= DomainMyLast; i++)
    {
      no = DomainNodeIndex[i];

      DomainTreeNodeLen[i] = Nodes[no].len;
    }

  for(level = 1; level < (1 << PTask); level++)
    {
      sendTask = ThisTask;
      recvTask = ThisTask ^ level;

      if(recvTask < NTask)
	MPI_Sendrecv(&DomainTreeNodeLen[DomainStartList[sendTask]],
		     (DomainEndList[sendTask] - DomainStartList[sendTask] + 1) * sizeof(FLOAT),
		     MPI_BYTE, recvTask, TAG_NODELEN,
		     &DomainTreeNodeLen[DomainStartList[recvTask]],
		     (DomainEndList[recvTask] - DomainStartList[recvTask] + 1) * sizeof(FLOAT),
		     MPI_BYTE, recvTask, TAG_NODELEN, MPI_COMM_WORLD, &status);
    }

  /* Finally, we update the top-level tree. */
  force_update_node_len_toptree();
}


/*! This function recursively enlarges nodes such that they always contain
 *  all their daughter nodes and daughter particles.
 */
void force_update_node_len_local(void)
{
  int i, p, k, no;
  FLOAT dist, distmax;

  for(i = 0; i < NumPart; i++)
    {
      no = Father[i];

      for(k = 0, distmax = 0; k < 3; k++)
	{
	  dist = P[i].Pos[k] - Nodes[no].center[k];
	  if(dist < 0)
	    dist = -dist;
	  if(dist > distmax)
	    distmax = dist;
	}

      if(distmax + distmax > Nodes[no].len)
	{
	  Nodes[no].len = distmax + distmax;
	  p = Nodes[no].u.d.father;

	  while(p >= 0)
	    {
	      distmax = Nodes[p].center[0] - Nodes[no].center[0];
	      if(distmax < 0)
		distmax = -distmax;
	      distmax = distmax + distmax + Nodes[no].len;

	      if(0.999999 * distmax > Nodes[p].len)
		{
		  Nodes[p].len = distmax;
		  no = p;
		  p = Nodes[p].u.d.father;
		}
	      else
		break;
	    }
	}
    }
}


/*! This function recursively enlarges nodes of the top-level tree such
 *  that they always contain all their daughter nodes.
 */
void force_update_node_len_toptree(void)
{
  int i, no, p;
  FLOAT distmax;

  for(i = 0; i < NTopleaves; i++)
    if(i < DomainMyStart || i > DomainMyLast)
      {
	no = DomainNodeIndex[i];

	if(Nodes[no].len < DomainTreeNodeLen[i])
	  Nodes[no].len = DomainTreeNodeLen[i];

	p = Nodes[no].u.d.father;

	while(p >= 0)
	  {
	    distmax = Nodes[p].center[0] - Nodes[no].center[0];
	    if(distmax < 0)
	      distmax = -distmax;
	    distmax = distmax + distmax + Nodes[no].len;

	    if(0.999999 * distmax > Nodes[p].len)
	      {
		Nodes[p].len = distmax;
		no = p;
		p = Nodes[p].u.d.father;
	      }
	    else
	      break;
	  }
      }
}




/*! This function updates the hmax-values in tree nodes that hold SPH
 *  particles. These values are needed to find all neighbors in the
 *  hydro-force computation.  Since the Hsml-values are potentially changed
 *  in the SPH-denity computation, force_update_hmax() should be carried
 *  out just before the hydrodynamical SPH forces are computed, i.e. after
 *  density().
 */
void force_update_hmax(void)
{
  int i, no;
  MPI_Status status;
  int level, sendTask, recvTask;

  force_update_node_hmax_local();

  for(i = DomainMyStart; i <= DomainMyLast; i++)
    {
      no = DomainNodeIndex[i];

      DomainHmax[i] = Extnodes[no].hmax;
    }

  /* share the hmax-data of the pseudo-particles accross CPUs */

  for(level = 1; level < (1 << PTask); level++)
    {
      sendTask = ThisTask;
      recvTask = ThisTask ^ level;

      if(recvTask < NTask)
	MPI_Sendrecv(&DomainHmax[DomainStartList[sendTask]],
		     (DomainEndList[sendTask] - DomainStartList[sendTask] + 1) * sizeof(FLOAT),
		     MPI_BYTE, recvTask, TAG_HMAX,
		     &DomainHmax[DomainStartList[recvTask]],
		     (DomainEndList[recvTask] - DomainStartList[recvTask] + 1) * sizeof(FLOAT),
		     MPI_BYTE, recvTask, TAG_HMAX, MPI_COMM_WORLD, &status);
    }


  force_update_node_hmax_toptree();
}

/*! This routine updates the hmax-values of local tree nodes.
 */
void force_update_node_hmax_local(void)
{
  int i, p, no;

  for(i = 0; i < N_gas; i++)
    {

      no = Father[i];

      if(SphP[i].Hsml > Extnodes[no].hmax)
	{

	  Extnodes[no].hmax = SphP[i].Hsml;
	  p = Nodes[no].u.d.father;

	  while(p >= 0)
	    {
	      if(Extnodes[no].hmax > Extnodes[p].hmax)
		{
		  Extnodes[p].hmax = Extnodes[no].hmax;
		  no = p;
		  p = Nodes[p].u.d.father;
		}
	      else
		break;
	    }
	}

    }
}




/*! This function recursively sets the hmax-values of the top-level tree.
 */
void force_update_node_hmax_toptree(void)
{

  int i, no, p;


  for(i = 0; i < NTopleaves; i++)
    if(i < DomainMyStart || i > DomainMyLast)
      {
	no = DomainNodeIndex[i];

	if(Extnodes[no].hmax < DomainHmax[i])
	  Extnodes[no].hmax = DomainHmax[i];

	p = Nodes[no].u.d.father;

	while(p >= 0)
	  {
	    if(Extnodes[no].hmax > Extnodes[p].hmax)
	      {
		Extnodes[p].hmax = Extnodes[no].hmax;
		no = p;
		p = Nodes[p].u.d.father;
	      }
	    else
	      break;
	  }
      }
}



/*! This routine computes the gravitational force for a given local
 *  particle, or for a particle in the communication buffer. Depending on
 *  the value of TypeOfOpeningCriterion, either the geometrical BH
 *  cell-opening criterion, or the `relative' opening criterion is used.
 */
int force_treeevaluate(int target, int mode, double *latticecountsum)
{
  struct NODE *nop = 0;
  int no, ninteractions, ptype;
  double r2[N_GRAVS], dx[N_GRAVS], dy[N_GRAVS], dz[N_GRAVS], mass[N_GRAVS], r, fac, h;
  double acc_x, acc_y, acc_z, pos_x, pos_y, pos_z, aold;
  
  int sG;
  FLOAT pmass;
  int i;
  double summass;
  double r2min, r2max;
  int pgravtype;

#ifdef NGRAVS_ACCUMULATOR
  long Nparticles[N_GRAVS];
#endif

#if defined(UNEQUALSOFTENINGS) && !defined(ADAPTIVE_GRAVSOFT_FORGAS)
  int maxsofttype;
#endif
#ifdef ADAPTIVE_GRAVSOFT_FORGAS
  double soft = 0;
#endif
#ifdef PERIODIC
  double boxsize, boxhalf;

  boxsize = All.BoxSize;
  boxhalf = 0.5 * All.BoxSize;
#endif

  acc_x = 0;
  acc_y = 0;
  acc_z = 0;
  ninteractions = 0;

  /* KC 8/9/14 mode seems to indicate whether the particle is local (0) or in the communication 
     buffer (1).  This needs to be extended to provide the mass of the target as well
     as, in general, we will now need this quantity
  */
  if(mode == 0)
    {
      pos_x = P[target].Pos[0];
      pos_y = P[target].Pos[1];
      pos_z = P[target].Pos[2];
      ptype = P[target].Type;
      pmass = P[target].Mass;
      aold = All.ErrTolForceAcc * P[target].OldAcc;
#ifdef ADAPTIVE_GRAVSOFT_FORGAS
      if(ptype == 0)
	soft = SphP[target].Hsml;
#endif
    }
  else
    {
      pos_x = GravDataGet[target].u.Pos[0];
      pos_y = GravDataGet[target].u.Pos[1];
      pos_z = GravDataGet[target].u.Pos[2];
    
      // KC 8/9/14 We always need this information now
      // specific to the particle under consideration!
      pmass = GravDataGet[target].Mass;
      ptype = GravDataGet[target].Type;

      aold = All.ErrTolForceAcc * GravDataGet[target].w.OldAcc;
#ifdef ADAPTIVE_GRAVSOFT_FORGAS
      if(ptype == 0)
	soft = GravDataGet[target].Soft;
#endif
    }



#ifndef UNEQUALSOFTENINGS
  h = All.ForceSoftening[ptype];
#endif
  no = All.MaxPart;		/* root node */

  while(no >= 0)
    {
      if(no < All.MaxPart)	/* single particle */
	{
	  /* the index of the node is the index of the particle */
	  /* observe the sign */
	  sG = TypeToGrav[P[no].Type];

	  mass[sG] = P[no].Mass;
	  dx[sG] = P[no].Pos[0] - pos_x;
	  dy[sG] = P[no].Pos[1] - pos_y;
	  dz[sG] = P[no].Pos[2] - pos_z;

#ifdef PERIODIC
	  dx[sG] = NEAREST(dx[sG]);
	  dy[sG] = NEAREST(dy[sG]);
	  dz[sG] = NEAREST(dz[sG]);
#endif
	  r2[sG] = dx[sG] * dx[sG] + dy[sG] * dy[sG] + dz[sG] * dz[sG];

#ifdef NGRAVS_ACCUMULATOR
	  Nparticles[sG] = 1;
#endif
      	}
      else
	{
	  if(no >= All.MaxPart + MaxNodes)	/* pseudo particle */
	    {
	      if(mode == 0)
		{
		  Exportflag[DomainTask[no - (All.MaxPart + MaxNodes)]] = 1;
		}
	      no = Nextnode[no - MaxNodes];
	      continue;
	    }
	  
	  r2min = INFINITY;
	  r2max = -INFINITY;
	  summass = 0;

	  nop = &Nodes[no];

	  for(i = 0; i < N_GRAVS; ++i) {
	    
#ifdef NGRAVS_ACCUMULATOR
	    Nparticles[i] = nop->u.d.Nparticles[i];
#endif
	    mass[i] = nop->u.d.mass[i];
	    summass += nop->u.d.mass[i];
	    dx[i] = nop->u.d.s[0][i] - pos_x;
	    dy[i] = nop->u.d.s[1][i] - pos_y;
	    dz[i] = nop->u.d.s[2][i] - pos_z;
	    
#ifdef PERIODIC
	    dx[i] = NEAREST(dx[i]);
	    dy[i] = NEAREST(dy[i]);
	    dz[i] = NEAREST(dz[i]);
#endif
	    r2[i] = dx[i] * dx[i] + dy[i] * dy[i] + dz[i] * dz[i];
	    
	    if(r2[i] < r2min)
	      r2min = r2[i];
	    if(r2[i] > r2max)
	      r2max = r2[i];
	    
	  }
	  
	  // KC 8/10/14 Flag that we will determine source type by the iteration over
	  // contribs
	  sG = -1;
	}

      // KC 8/14/14 This is the single particle criterion check...
      if(no < All.MaxPart)
	{
#ifdef UNEQUALSOFTENINGS
#ifdef ADAPTIVE_GRAVSOFT_FORGAS
	  if(ptype == 0)
	    h = soft;
	  else
	    h = All.ForceSoftening[ptype];

	  if(P[no].Type == 0)
	    {
	      if(h < SphP[no].Hsml)
		h = SphP[no].Hsml;
	    }
	  else
	    {
	      if(h < All.ForceSoftening[P[no].Type])
		h = All.ForceSoftening[P[no].Type];
	    }
#else
	  h = All.ForceSoftening[ptype];
	  if(h < All.ForceSoftening[P[no].Type])
	    h = All.ForceSoftening[P[no].Type];
#endif
#endif
	  no = Nextnode[no];
	}
      else			/* we have an  internal node. Need to check opening criterion */
	{
	  if(mode == 1)
	    {
	      if((nop->u.d.bitflags & 3) == 1)	/* if it's a top-level node
						 * which does not contain
						 * local particles we can
						 * continue to do a short-cut */
		{
		  no = nop->u.d.sibling;
		  continue;
		}
	    }


	  if(All.ErrTolTheta)	/* check Barnes-Hut opening criterion */
	    {
	      if(nop->len * nop->len > r2min * All.ErrTolTheta * All.ErrTolTheta)
		{
		  /* open cell */
		  no = nop->u.d.nextnode;
		  continue;
		}
	    }
	  else			/* check relative opening criterion */
	    {
	      // KC 11/3/14
	      // The computation here avoids divisions.  It is:
	      //
	      // M/r^2 (l/r)^2 > aold
	      //
	      // All interactions are bounded above by Newton, so this check is fine
	      if(summass * nop->len * nop->len > r2min * r2min * aold)
		{
		  /* open cell */
		  no = nop->u.d.nextnode;
		  continue;
		}

	      /* check in addition whether we lie inside the cell */
	      if(fabs(nop->center[0] - pos_x) < 0.60 * nop->len)
		{
		  if(fabs(nop->center[1] - pos_y) < 0.60 * nop->len)
		    {
		      if(fabs(nop->center[2] - pos_z) < 0.60 * nop->len)
			{
			  no = nop->u.d.nextnode;
			  continue;
			}
		    }
		}
	    }

#ifdef UNEQUALSOFTENINGS
#ifndef ADAPTIVE_GRAVSOFT_FORGAS
	  h = All.ForceSoftening[ptype];
	  maxsofttype = (nop->u.d.bitflags >> 2) & 7;
	  if(maxsofttype == 7) /* may only occur for zero mass top-level nodes */
	    {
	      if(summass > 0)
		endrun(986);
	      no = nop->u.d.nextnode;
	      continue;
	    }
	  else
	    {
	      if(h < All.ForceSoftening[maxsofttype])
		{
		  h = All.ForceSoftening[maxsofttype];
		  
		  if(r2max < h * h)
		    {
		      if(((nop->u.d.bitflags >> 5) & 1))	/* bit-5 signals that there are particles of different softening in the node */
			{
			  no = nop->u.d.nextnode;
			  continue;
			}
		    }
		}
	    }
#else
	  if(ptype == 0)
	    h = soft;
	  else
	    h = All.ForceSoftening[ptype];

	  if(h < nop->maxsoft)
	    {
	      h = nop->maxsoft;
	      if(r2max < h * h)
		{
		  no = nop->u.d.nextnode;
		  continue;
		}
	    }
#endif
#endif

	  no = nop->u.d.sibling;	/* ok, node can be used */

	  if(mode == 1)
	    {
	      if(((nop->u.d.bitflags) & 1))	/* Bit 0 signals that this node belongs to top-level tree */
		continue;
	    }
	}
    
      /* KC 8/9/14 Here we sum over all the contributions from the various gravitational 
	 constituents of the node. If its a 
	 pseudo-particle, we have to sum over all the different contributions
      */

      pgravtype = TypeToGrav[ptype];

      if(sG > -1) {
	
	// A single particle contributes
	r = sqrt(r2[sG]);
	  
	// KC 1/31/15
	// AccelFxns now return the normalized value in one step...
	if(r >= h)
	  fac = (*AccelFxns[pgravtype][sG])(pmass, mass[sG], h, r, 1);
	else
	  fac = (*AccelSplines[pgravtype][sG])(pmass, mass[sG], h, r, 1);
	
	acc_x += dx[sG] * fac;
	acc_y += dy[sG] * fac;
	acc_z += dz[sG] * fac;
      }
      else {

	// A node contributes
	for(i = 0; i < N_GRAVS; ++i) {

	  if(mass[i] == 0.0)
	    continue;

	  // KC 8/13/14 This code now uses i appropriately
	  r = sqrt(r2[i]);
	  
	  // KC 10/18/14
	  if(r >= h) {
#ifdef NGRAVS_ACCUMULATOR
	    fac = (*AccelFxns[pgravtype][i])(pmass, mass[i], h, r, Nparticles[i]);
#else
	    fac = (*AccelFxns[pgravtype][i])(pmass, mass[i], h, r, 1);
#endif
	  }
	  else
	    {
#ifdef NGRAVS_ACCUMULATOR
	      fac = (*AccelSplines[pgravtype][i])(pmass, mass[i], h, r, Nparticles[i]);
#else
	      fac = (*AccelSplines[pgravtype][i])(pmass, mass[i], h, r, 1);
#endif 
	    }
	  
	  acc_x += dx[i] * fac;
	  acc_y += dy[i] * fac;
	  acc_z += dz[i] * fac;
	  
	} /* KC 8/9/14 closes sum over contributions */
      }
      
      ninteractions++;
    }


  /* store result at the proper place */
  if(mode == 0)
    {
      P[target].GravAccel[0] = acc_x;
      P[target].GravAccel[1] = acc_y;
      P[target].GravAccel[2] = acc_z;
      P[target].GravCost = ninteractions;
    }
  else
    {
      GravDataResult[target].u.Acc[0] = acc_x;
      GravDataResult[target].u.Acc[1] = acc_y;
      GravDataResult[target].u.Acc[2] = acc_z;
      GravDataResult[target].w.Ninteractions = ninteractions;
    }

#ifdef PERIODIC
  *latticecountsum += force_treeevaluate_lattice_correction(target, mode, pos_x, pos_y, pos_z, aold);
#endif

  return ninteractions;
}

#ifdef PMGRID
/*! In the TreePM algorithm, the tree is walked only locally around the
 *  target coordinate.  Tree nodes that fall outside a box of half
 *  side-length Rcut= RCUT*ASMTH*MeshSize can be discarded. The short-range
 *  potential is modified by a complementary error function, multiplied
 *  with the Newtonian form. The resulting short-range suppression compared
 *  to the Newtonian force is tabulated, because looking up from this table
 *  is faster than recomputing the corresponding factor, despite the
 *  memory-access panelty (which reduces cache performance) incurred by the
 *  table.
 */
int force_treeevaluate_shortrange(int target, int mode)
{

  struct NODE *nop = 0;
  int no, ptype, ninteractions, tabindex;
  double r2[N_GRAVS], dx[N_GRAVS], dy[N_GRAVS], dz[N_GRAVS], mass[N_GRAVS], r, fac, h;
  double acc_x, acc_y, acc_z, pos_x, pos_y, pos_z, aold;
  double eff_dist;
  double rcut, asmth, asmthfac, rcut2, dist;

  // KC 10/17/14
  // ngravs extensions
  int sG;
  FLOAT pmass;
  int i;
  double summass;
  double r2min, r2max;
  char nintflag;
  int pgravtype, whichGrav;

#ifdef NGRAVS_ACCUMULATOR
  long Nparticles[N_GRAVS];
#endif

#if defined(UNEQUALSOFTENINGS) && !defined(ADAPTIVE_GRAVSOFT_FORGAS)
  int maxsofttype;
#endif
#ifdef ADAPTIVE_GRAVSOFT_FORGAS
  double soft = 0;
#endif
#ifdef PERIODIC
  double boxsize, boxhalf;

  boxsize = All.BoxSize;
  boxhalf = 0.5 * All.BoxSize;
#endif

  acc_x = 0;
  acc_y = 0;
  acc_z = 0;
  ninteractions = 0;

  // KC 10/17/14
  // Is particle local or in the comm buffer?
  if(mode == 0)
    {
      pos_x = P[target].Pos[0];
      pos_y = P[target].Pos[1];
      pos_z = P[target].Pos[2];
      ptype = P[target].Type;
      pmass = P[target].Mass;
      aold = All.ErrTolForceAcc * P[target].OldAcc;
#ifdef ADAPTIVE_GRAVSOFT_FORGAS
      if(ptype == 0)
	soft = SphP[target].Hsml;
#endif
    }
  else
    {
      pos_x = GravDataGet[target].u.Pos[0];
      pos_y = GravDataGet[target].u.Pos[1];
      pos_z = GravDataGet[target].u.Pos[2];

      ptype = GravDataGet[target].Type;
      pmass = GravDataGet[target].Mass;

      aold = All.ErrTolForceAcc * GravDataGet[target].w.OldAcc;
#ifdef ADAPTIVE_GRAVSOFT_FORGAS
      if(ptype == 0)
	soft = GravDataGet[target].Soft;
#endif
    }

  rcut = All.Rcut[0];
  asmth = All.Asmth[0];
#ifdef PLACEHIGHRESREGION
  if(((1 << ptype) & (PLACEHIGHRESREGION)))
    {
      rcut = All.Rcut[1];
      asmth = All.Asmth[1];
    }
#endif
  rcut2 = rcut * rcut;

  asmthfac = 0.5 / asmth * (NTAB / 3.0);

  // KC 1/31/15
  // We will be computing these in the individual splines anyway
#ifndef UNEQUALSOFTENINGS
  h = All.ForceSoftening[ptype];
#endif

  no = All.MaxPart;		/* root node */

  while(no >= 0)
    {
      if(no < All.MaxPart)
	{

	  sG = TypeToGrav[P[no].Type];

	  /* the index of the node is the index of the particle */
	  mass[sG] = P[no].Mass;
	  dx[sG] = P[no].Pos[0] - pos_x;
	  dy[sG] = P[no].Pos[1] - pos_y;
	  dz[sG] = P[no].Pos[2] - pos_z;
#ifdef PERIODIC
	  dx[sG] = NEAREST(dx[sG]);
	  dy[sG] = NEAREST(dy[sG]);
	  dz[sG] = NEAREST(dz[sG]);
#endif
	  r2[sG] = dx[sG] * dx[sG] + dy[sG] * dy[sG] + dz[sG] * dz[sG];

#ifdef NGRAVS_ACCUMULATOR
	  Nparticles[sG] = 1;
#endif
	
#ifdef UNEQUALSOFTENINGS
#ifdef ADAPTIVE_GRAVSOFT_FORGAS
	  if(ptype == 0)
	    h = soft;
	  else
	    h = All.ForceSoftening[ptype];

	  if(P[no].Type == 0)
	    {
	      if(h < SphP[no].Hsml)
		h = SphP[no].Hsml;
	    }
	  else
	    {
	      if(h < All.ForceSoftening[P[no].Type])
		h = All.ForceSoftening[P[no].Type];
	    }
#else
	  h = All.ForceSoftening[ptype];
	  if(h < All.ForceSoftening[P[no].Type])
	    h = All.ForceSoftening[P[no].Type];
#endif
#endif
	  no = Nextnode[no];
	}
      else			/* we have an  internal node */
	{
	  if(no >= All.MaxPart + MaxNodes)	/* pseudo particle */
	    {
	      if(mode == 0)
		{
		  Exportflag[DomainTask[no - (All.MaxPart + MaxNodes)]] = 1;
		}
	      no = Nextnode[no - MaxNodes];
	      continue;
	    }
	  
	  nop = &Nodes[no];
	  
	  if(mode == 1)
	    {
	      if((nop->u.d.bitflags & 3) == 1)	/* if it's a top-level node
						 * which does not contain
						 * local particles we can
						 * continue at this point
						 */
		{
		  no = nop->u.d.sibling;
		  continue;
		}
	    }

	  /* KC 10/17/14 A pseudo-particle now represents N_GRAV distinct
	     pseudo-particles, so we need the arrays */
	  r2min = INFINITY;
	  r2max = -INFINITY;
	  summass = 0;
	  
	  for(i = 0; i < N_GRAVS; ++i) {
	  
#ifdef NGRAVS_ACCUMULATOR
	    Nparticles[i] = nop->u.d.Nparticles[i];
#endif
	    mass[i] = nop->u.d.mass[i];
	    summass += nop->u.d.mass[i];
	    dx[i] = nop->u.d.s[0][i] - pos_x;
	    dy[i] = nop->u.d.s[1][i] - pos_y;
	    dz[i] = nop->u.d.s[2][i] - pos_z;

#ifdef PERIODIC
	    dx[i] = NEAREST(dx[i]);
	    dy[i] = NEAREST(dy[i]);
	    dz[i] = NEAREST(dz[i]);
#endif
	    r2[i] = dx[i] * dx[i] + dy[i] * dy[i] + dz[i] * dz[i];
	    
	    if(r2[i] < r2min)
	      r2min = r2[i];
	    if(r2[i] > r2max)
	      r2max = r2[i];
	  }
	
	  sG = -1;

	  if(r2min > rcut2)
	    {
	      /* check whether we can stop walking along this branch */
	      eff_dist = rcut + 0.5 * nop->len;
#ifdef PERIODIC
	      dist = NEAREST(nop->center[0] - pos_x);
#else
	      dist = nop->center[0] - pos_x;
#endif
	      if(dist < -eff_dist || dist > eff_dist)
		{
		  no = nop->u.d.sibling;
		  continue;
		}
#ifdef PERIODIC
	      dist = NEAREST(nop->center[1] - pos_y);
#else
	      dist = nop->center[1] - pos_y;
#endif
	      if(dist < -eff_dist || dist > eff_dist)
		{
		  no = nop->u.d.sibling;
		  continue;
		}
#ifdef PERIODIC
	      dist = NEAREST(nop->center[2] - pos_z);
#else
	      dist = nop->center[2] - pos_z;
#endif
	      if(dist < -eff_dist || dist > eff_dist)
		{
		  no = nop->u.d.sibling;
		  continue;
		}
	    }


	  if(All.ErrTolTheta)	/* check Barnes-Hut opening criterion */
	    {
	      if(nop->len * nop->len > r2min * All.ErrTolTheta * All.ErrTolTheta)
		{
		  /* open cell */
		  no = nop->u.d.nextnode;
		  continue;
		}
	    }
	  else			/* check relative opening criterion */
	    {
	      
	      if(summass * nop->len * nop->len > r2min * r2min * aold)
		{
		  /* open cell */
		  no = nop->u.d.nextnode;
		  continue;
		}

	      /* check in addition whether we lie inside the cell */

	      if(fabs(nop->center[0] - pos_x) < 0.60 * nop->len)
		{
		  if(fabs(nop->center[1] - pos_y) < 0.60 * nop->len)
		    {
		      if(fabs(nop->center[2] - pos_z) < 0.60 * nop->len)
			{
			  no = nop->u.d.nextnode;
			  continue;
			}
		    }
		}
	    }

#ifdef UNEQUALSOFTENINGS
#ifndef ADAPTIVE_GRAVSOFT_FORGAS
	  h = All.ForceSoftening[ptype];
          maxsofttype = (nop->u.d.bitflags >> 2) & 7;
          if(maxsofttype == 7) /* may only occur for zero mass top-level nodes */
            {
              
	      if(summass > 0)
                endrun(987);
              no = nop->u.d.nextnode;
              continue;
            }
          else
            {
              if(h < All.ForceSoftening[maxsofttype])
                {
                  h = All.ForceSoftening[maxsofttype];
                  if(r2max < h * h)
                    {
                      if(((nop->u.d.bitflags >> 5) & 1))	/* bit-5 signals that there are particles of different softening in the node */
                        {
                          no = nop->u.d.nextnode;
                          
                          continue;
                        }
                    }
                }
            }
#else
	  if(ptype == 0)
	    h = soft;
	  else
	    h = All.ForceSoftening[ptype];

	  if(h < nop->maxsoft)
	    {
	      h = nop->maxsoft;
	      if(r2max < h * h)
		{
		  no = nop->u.d.nextnode;
		  continue;
		}
	    }
#endif
#endif
	  no = nop->u.d.sibling;	/* ok, node can be used */

	  if(mode == 1)
	    {
	      if(((nop->u.d.bitflags) & 1))	/* Bit 0 signals that this node belongs to top-level tree */
		continue;
	    }
	}

      // KC 10/17/14
      // Now iterate over the N_GRAVS and do what's right
      nintflag = 0;
      pgravtype = TypeToGrav[ptype];

      if(sG > -1) {

	 r = sqrt(r2[sG]);

	 if(r >= h)
	   fac = (*AccelFxns[pgravtype][sG])(pmass, mass[sG], h, r, 1);
	 else 
	   fac = (*AccelSplines[pgravtype][sG])(pmass, mass[sG], h, r, 1);
	 
	 tabindex = (int) (asmthfac * r);

	  if(tabindex < NTAB)
	    {
	      fac *= shortrange_table[tabindex];

	      acc_x += dx[sG] * fac;
	      acc_y += dy[sG] * fac;
	      acc_z += dz[sG] * fac;
	      
	      // Flag to record interactions
	      nintflag = 1;
	    }
      }
      else {
	
	for(whichGrav = 0; whichGrav < N_GRAVS; ++whichGrav) {

	  if(mass[whichGrav] == 0.0)
	    continue;

	  r = sqrt(r2[whichGrav]);

	  if(r >= h) {
#ifdef NGRAVS_ACCUMULATOR

	    fac = (*AccelFxns[pgravtype][whichGrav])(pmass, mass[whichGrav], h, r, Nparticles[whichGrav]);
#else
	    fac = (*AccelFxns[pgravtype][whichGrav])(pmass, mass[whichGrav], h, r, 1);
#endif
	  }
	  else {
#ifdef NGRAVS_ACCUMULATOR
	    fac = (*AccelSplines[pgravtype][whichGrav])(pmass, mass[whichGrav], h, r, Nparticles[whichGrav]);
#else
	    fac = (*AccelSplines[pgravtype][whichGrav])(pmass, mass[whichGrav], h, r, 1);
#endif
	  }
	  tabindex = (int) (asmthfac * r);

	  if(tabindex < NTAB)
	    {
	      fac *= shortrange_table[tabindex];

	      acc_x += dx[whichGrav] * fac;
	      acc_y += dy[whichGrav] * fac;
	      acc_z += dz[whichGrav] * fac;
	      
	      // Flag to record interactions
	      nintflag = 1;
	    }
	}      
      }

      // KC 10/17/14
      // For consistency with extension of force_treeevaluate()
      if(nintflag)
	ninteractions++;
    }

  /* store result at the proper place */
  if(mode == 0)
    {
      P[target].GravAccel[0] = acc_x;
      P[target].GravAccel[1] = acc_y;
      P[target].GravAccel[2] = acc_z;
      P[target].GravCost = ninteractions;
    }
  else
    {
      GravDataResult[target].u.Acc[0] = acc_x;
      GravDataResult[target].u.Acc[1] = acc_y;
      GravDataResult[target].u.Acc[2] = acc_z;
      GravDataResult[target].w.Ninteractions = ninteractions;
    }

  return ninteractions;
}

#endif



#if defined PERIODIC
/*! This function computes the Lattice sum correction, and is needed if periodic
 *  boundary conditions together with a pure tree algorithm are used. Note
 *  that the ordinary tree walk does not carry out this correction directly
 *  as it was done in Gadget-1.1. Instead, the tree is walked a second
 *  time. This is actually faster because the "Lattice sum-Treewalk" can use a
 *  different opening criterion than the normal tree walk. In particular,
 *  the Lattice sum correction is negligible for particles that are very close,
 *  but it is large for particles that are far away (this is quite
 *  different for the normal direct force). So we can here use a different
 *  opening criterion. Sufficient accuracy is usually obtained if the node
 *  length has dropped to a certain fraction ~< 0.25 of the
 *  BoxLength. However, we may only short-cut the interaction list of the
 *  normal full Lattice sum tree walk if we are sure that the whole node and all
 *  daughter nodes "lie on the same side" of the periodic boundary,
 *  i.e. that the real tree walk would not find a daughter node or particle
 *  that was mapped to a different nearest neighbour position when the tree
 *  walk would be further refined.
 */
int force_treeevaluate_lattice_correction(int target, int mode, double pos_x, double pos_y, double pos_z,
					double aold)
{
  struct NODE *nop = 0;
  int no, cost;
  double dx[N_GRAVS], dy[N_GRAVS], dz[N_GRAVS], mass[N_GRAVS], r2[N_GRAVS];
  int signx, signy, signz;
  int i, j, k, openflag;
  double u, v, w;
  double f1, f2, f3, f4, f5, f6, f7, f8;
  double acc_x, acc_y, acc_z;
  double boxsize, boxhalf;

  // KC 10/25/14
  // N_GRAVS variables
  int sG;
  int n;
  double summass;
  double r2min, r2max;
  int pgravtype;

  boxsize = All.BoxSize;
  boxhalf = 0.5 * All.BoxSize;

  acc_x = 0;
  acc_y = 0;
  acc_z = 0;
  cost = 0;

  no = All.MaxPart;

  while(no >= 0)
    {
      if(no < All.MaxPart)	/* single particle */
	{
	  /* the index of the node is the index of the particle */
	  /* observe the sign */

	  sG = TypeToGrav[P[no].Type];
	
	  mass[sG] = P[no].Mass;
	  dx[sG] = NEAREST(P[no].Pos[0] - pos_x);
	  dy[sG] = NEAREST(P[no].Pos[1] - pos_y);
	  dz[sG] = NEAREST(P[no].Pos[2] - pos_z);
	  
	  r2[sG] = dx[sG] * dx[sG] + dy[sG] * dy[sG] + dz[sG] * dz[sG];
	}
      else
	{
	  if(no >= All.MaxPart + MaxNodes)	/* pseudo particle */
	    {
	      if(mode == 0)
		{
		  Exportflag[DomainTask[no - (All.MaxPart + MaxNodes)]] = 1;
		}

	      no = Nextnode[no - MaxNodes];
	      continue;
	    }

	  nop = &Nodes[no];

	  r2min = INFINITY;
	  r2max = -INFINITY;
	  summass = 0;

	  for(i = 0; i < N_GRAVS; ++i) {
	    
	    mass[i] = nop->u.d.mass[i];
	    summass += nop->u.d.mass[i];
	    dx[i] = NEAREST(nop->u.d.s[0][i] - pos_x);
	    dy[i] = NEAREST(nop->u.d.s[1][i] - pos_y);
	    dz[i] = NEAREST(nop->u.d.s[2][i] - pos_z);
	    r2[i] = dx[i] * dx[i] + dy[i] * dy[i] + dz[i] * dz[i];
	    
	    if(r2[i] < r2min)
	      r2min = r2[i];
	    if(r2[i] > r2max)
	      r2max = r2[i];
	  }
	  
	  sG = -1;
	}

      if(no < All.MaxPart)
	no = Nextnode[no];
      else			/* we have an  internal node. Need to check opening criterion */
	{
	  openflag = 0;

	  // KC 10/25/14
	  // SSS We don't need r2 unless we are opening, but we have already computed it anyway...
	  // Again, crudely go minimum
	  if(All.ErrTolTheta)	/* check Barnes-Hut opening criterion */
	    {
	      if(nop->len * nop->len > r2min * All.ErrTolTheta * All.ErrTolTheta)
		{
		  openflag = 1;
		}
	    }
	  else			/* check relative opening criterion */
	    {
	      if(summass * nop->len * nop->len > r2min * r2min * aold)
		{
		  openflag = 1;
		}
	      else
		{
		  if(fabs(nop->center[0] - pos_x) < 0.60 * nop->len)
		    {
		      if(fabs(nop->center[1] - pos_y) < 0.60 * nop->len)
			{
			  if(fabs(nop->center[2] - pos_z) < 0.60 * nop->len)
			    {
			      openflag = 1;
			    }
			}
		    }
		}
	    }

	  if(openflag)
	    {
	      /* now we check if we can avoid opening the cell */

	      u = nop->center[0] - pos_x;
	      if(u > boxhalf)
		u -= boxsize;
	      if(u < -boxhalf)
		u += boxsize;

	      if(fabs(u) > 0.5 * (boxsize - nop->len))
		{
		  no = nop->u.d.nextnode;
		  continue;
		}

	      u = nop->center[1] - pos_y;
	      if(u > boxhalf)
		u -= boxsize;
	      if(u < -boxhalf)
		u += boxsize;

	      if(fabs(u) > 0.5 * (boxsize - nop->len))
		{
		  no = nop->u.d.nextnode;
		  continue;
		}

	      u = nop->center[2] - pos_z;
	      if(u > boxhalf)
		u -= boxsize;
	      if(u < -boxhalf)
		u += boxsize;

	      if(fabs(u) > 0.5 * (boxsize - nop->len))
		{
		  no = nop->u.d.nextnode;
		  continue;
		}

	      /* if the cell is too large, we need to refine
	       * it further 
	       */
	      if(nop->len > 0.20 * boxsize)
		{
		  /* cell is too large */
		  no = nop->u.d.nextnode;
		  continue;
		}
	    }

	  no = nop->u.d.sibling;	/* ok, node can be used */

	  if(mode == 1)
	    {
	      if((nop->u.d.bitflags & 1))	/* Bit 0 signals that this node belongs to top-level tree */
		continue;
	    }
	}

      // KC 12/6/14
      // Compute the lattice correction force due to everything everywhere
      pgravtype = TypeToGrav[P[target].Type];

      if(sG < 0) {
	
	// KC 12/6/14
	// Iterate over all contributions from the node
 	for(n = 0; n < N_GRAVS; ++n) {

	  if(mass[n] == 0.0)
	    continue;

	  if(dx[n] < 0)
	    {
	      dx[n] = -dx[n];
	      signx = +1;
	    }
	  else
	    signx = -1;

	  if(dy[n] < 0)
	    {
	      dy[n] = -dy[n];
	      signy = +1;
	    }
	  else
	    signy = -1;

	  if(dz[n] < 0)
	    {
	      dz[n] = -dz[n];
	      signz = +1;
	    }
	  else
	    signz = -1;

	  u = dx[n] * fac_intp;
	  i = (int) u;
	  if(i >= EN)
	    i = EN - 1;
	  u -= i;
	  v = dy[n] * fac_intp;
	  j = (int) v;
	  if(j >= EN)
	    j = EN - 1;
	  v -= j;
	  w = dz[n] * fac_intp;
	  k = (int) w;
	  if(k >= EN)
	    k = EN - 1;
	  w -= k;

	  /* compute factors for trilinear interpolation */

	  f1 = (1 - u) * (1 - v) * (1 - w);
	  f2 = (1 - u) * (1 - v) * (w);
	  f3 = (1 - u) * (v) * (1 - w);
	  f4 = (1 - u) * (v) * (w);
	  f5 = (u) * (1 - v) * (1 - w);
	  f6 = (u) * (1 - v) * (w);
	  f7 = (u) * (v) * (1 - w);
	  f8 = (u) * (v) * (w);

	  // KC 12/6/14
	  // Note that the force correction computed and stored in fcorr is for a point source.
	  // The way that a node's contribution is reckoned is by taking this point contribution, 
	  // multiplying it by the mass, but dividing out by the number of contributors within that
	  // node when determining scale.  Here, the scale applied is always that of a single 
	  // contributor, WITH NO MASS PREFACTOR.
	  acc_x += mass[n] * signx * (fcorrx[pgravtype][n][i][j][k] * f1 +
				      fcorrx[pgravtype][n][i][j][k + 1] * f2 +
				      fcorrx[pgravtype][n][i][j + 1][k] * f3 +
				      fcorrx[pgravtype][n][i][j + 1][k + 1] * f4 +
				      fcorrx[pgravtype][n][i + 1][j][k] * f5 +
				      fcorrx[pgravtype][n][i + 1][j][k + 1] * f6 +
				      fcorrx[pgravtype][n][i + 1][j + 1][k] * f7 + fcorrx[pgravtype][n][i + 1][j + 1][k + 1] * f8);

	  acc_y += mass[n] * signy * (fcorry[pgravtype][n][i][j][k] * f1 +
				      fcorry[pgravtype][n][i][j][k + 1] * f2 +
				      fcorry[pgravtype][n][i][j + 1][k] * f3 +
				      fcorry[pgravtype][n][i][j + 1][k + 1] * f4 +
				      fcorry[pgravtype][n][i + 1][j][k] * f5 +
				      fcorry[pgravtype][n][i + 1][j][k + 1] * f6 +
				      fcorry[pgravtype][n][i + 1][j + 1][k] * f7 + fcorry[pgravtype][n][i + 1][j + 1][k + 1] * f8);

	  acc_z += mass[n] * signz * (fcorrz[pgravtype][n][i][j][k] * f1 +
				      fcorrz[pgravtype][n][i][j][k + 1] * f2 +
				      fcorrz[pgravtype][n][i][j + 1][k] * f3 +
				      fcorrz[pgravtype][n][i][j + 1][k + 1] * f4 +
				      fcorrz[pgravtype][n][i + 1][j][k] * f5 +
				      fcorrz[pgravtype][n][i + 1][j][k + 1] * f6 +
				      fcorrz[pgravtype][n][i + 1][j + 1][k] * f7 + fcorrz[pgravtype][n][i + 1][j + 1][k + 1] * f8);
	}
      }
      else {

	if(dx[sG] < 0)
	  {
	    dx[sG] = -dx[sG];
	    signx = +1;
	  }
	else
	  signx = -1;

	if(dy[sG] < 0)
	  {
	    dy[sG] = -dy[sG];
	    signy = +1;
	  }
	else
	  signy = -1;

	if(dz[sG] < 0)
	  {
	    dz[sG] = -dz[sG];
	    signz = +1;
	  }
	else
	  signz = -1;

	u = dx[sG] * fac_intp;
	i = (int) u;
	if(i >= EN)
	  i = EN - 1;
	u -= i;
	v = dy[sG] * fac_intp;
	j = (int) v;
	if(j >= EN)
	  j = EN - 1;
	v -= j;
	w = dz[sG] * fac_intp;
	k = (int) w;
	if(k >= EN)
	  k = EN - 1;
	w -= k;

	/* compute factors for trilinear interpolation */

	f1 = (1 - u) * (1 - v) * (1 - w);
	f2 = (1 - u) * (1 - v) * (w);
	f3 = (1 - u) * (v) * (1 - w);
	f4 = (1 - u) * (v) * (w);
	f5 = (u) * (1 - v) * (1 - w);
	f6 = (u) * (1 - v) * (w);
	f7 = (u) * (v) * (1 - w);
	f8 = (u) * (v) * (w);

	acc_x += mass[sG] * signx * (fcorrx[pgravtype][sG][i][j][k] * f1 +
				    fcorrx[pgravtype][sG][i][j][k + 1] * f2 +
				    fcorrx[pgravtype][sG][i][j + 1][k] * f3 +
				    fcorrx[pgravtype][sG][i][j + 1][k + 1] * f4 +
				    fcorrx[pgravtype][sG][i + 1][j][k] * f5 +
				    fcorrx[pgravtype][sG][i + 1][j][k + 1] * f6 +
				    fcorrx[pgravtype][sG][i + 1][j + 1][k] * f7 + fcorrx[pgravtype][sG][i + 1][j + 1][k + 1] * f8);

	acc_y += mass[sG] * signy * (fcorry[pgravtype][sG][i][j][k] * f1 +
				    fcorry[pgravtype][sG][i][j][k + 1] * f2 +
				    fcorry[pgravtype][sG][i][j + 1][k] * f3 +
				    fcorry[pgravtype][sG][i][j + 1][k + 1] * f4 +
				    fcorry[pgravtype][sG][i + 1][j][k] * f5 +
				    fcorry[pgravtype][sG][i + 1][j][k + 1] * f6 +
				    fcorry[pgravtype][sG][i + 1][j + 1][k] * f7 + fcorry[pgravtype][sG][i + 1][j + 1][k + 1] * f8);

	acc_z += mass[sG] * signz * (fcorrz[pgravtype][sG][i][j][k] * f1 +
				    fcorrz[pgravtype][sG][i][j][k + 1] * f2 +
				    fcorrz[pgravtype][sG][i][j + 1][k] * f3 +
				    fcorrz[pgravtype][sG][i][j + 1][k + 1] * f4 +
				    fcorrz[pgravtype][sG][i + 1][j][k] * f5 +
				    fcorrz[pgravtype][sG][i + 1][j][k + 1] * f6 +
				    fcorrz[pgravtype][sG][i + 1][j + 1][k] * f7 + fcorrz[pgravtype][sG][i + 1][j + 1][k + 1] * f8);
      }
      cost++;
    }


  /* add the result at the proper place */

  if(mode == 0)
    {
      P[target].GravAccel[0] += acc_x;
      P[target].GravAccel[1] += acc_y;
      P[target].GravAccel[2] += acc_z;
      P[target].GravCost += cost;
    }
  else
    {
      GravDataResult[target].u.Acc[0] += acc_x;
      GravDataResult[target].u.Acc[1] += acc_y;
      GravDataResult[target].u.Acc[2] += acc_z;
      GravDataResult[target].w.Ninteractions += cost;
    }

  return cost;
}

#endif





/*! This routine computes the gravitational potential by walking the
 *  tree. The same opening criteria is used as for the gravitational force
 *  walk.
 */
void force_treeevaluate_potential(int target, int mode)
{
#ifdef OUTPUTPOTENTIAL
  struct NODE *nop = 0;
  int no, ptype;
  double r2[N_GRAVS], dx[N_GRAVS], dy[N_GRAVS], dz[N_GRAVS], mass[N_GRAVS], r, h, h_inv;
  double pot, pos_x, pos_y, pos_z, aold;
#if defined(UNEQUALSOFTENINGS) && !defined(ADAPTIVE_GRAVSOFT_FORGAS)
  int maxsofttype;
#endif
#ifdef ADAPTIVE_GRAVSOFT_FORGAS
  double soft = 0;
#endif
#ifdef PERIODIC
  double boxsize, boxhalf;
  boxsize = All.BoxSize;
  boxhalf = 0.5 * All.BoxSize;
#endif

  // ngravs extension
  double summass, r2min, r2max;
  FLOAT pmass;
  int i, pgravtype, sG;
#ifdef NGRAVS_ACCUMULATOR
  long Nparticles[N_GRAVS];
#endif


  pot = 0;

  if(mode == 0)
    {
      pos_x = P[target].Pos[0];
      pos_y = P[target].Pos[1];
      pos_z = P[target].Pos[2];
      ptype = P[target].Type;
      pmass = P[target].Mass;

      aold = All.ErrTolForceAcc * P[target].OldAcc;
#ifdef ADAPTIVE_GRAVSOFT_FORGAS
      if(ptype == 0)
	soft = SphP[target].Hsml;
#endif
    }
  else
    {
      pos_x = GravDataGet[target].u.Pos[0];
      pos_y = GravDataGet[target].u.Pos[1];
      pos_z = GravDataGet[target].u.Pos[2];
      ptype = GravDataGet[target].Type;
      pmass = GravDataGet[target].Mass;

      aold = All.ErrTolForceAcc * GravDataGet[target].w.OldAcc;
#ifdef ADAPTIVE_GRAVSOFT_FORGAS
      if(ptype == 0)
	soft = GravDataGet[target].Soft;
#endif
    }


#ifndef UNEQUALSOFTENINGS
  h = All.ForceSoftening[ptype];
#endif
  no = All.MaxPart;

  while(no >= 0)
    {
      if(no < All.MaxPart)	/* single particle */
	{
	  /* the index of the node is the index of the particle */
	  /* observe the sign */
	  
	  sG = TypeToGrav[P[no].Type];
	  mass[sG] = P[no].Mass;
	  dx[sG] = P[no].Pos[0] - pos_x;
	  dy[sG] = P[no].Pos[1] - pos_y;
	  dz[sG] = P[no].Pos[2] - pos_z;

#ifdef PERIODIC
	  dx[sG] = NEAREST(dx[sG]);
	  dy[sG] = NEAREST(dy[sG]);
	  dz[sG] = NEAREST(dz[sG]);
#endif
	  r2[sG] = dx[sG] * dx[sG] + dy[sG] * dy[sG] + dz[sG] * dz[sG];
	}
      else
	{
	  if(no >= All.MaxPart + MaxNodes)	/* pseudo particle */
	    {
	      if(mode == 0)
		{
		  Exportflag[DomainTask[no - (All.MaxPart + MaxNodes)]] = 1;
		}
	      no = Nextnode[no - MaxNodes];
	      continue;
	    }

	  r2min = INFINITY;
	  r2max = -INFINITY;
	  summass = 0;

	  nop = &Nodes[no];

	  for(i = 0; i < N_GRAVS; ++i) {
#ifdef NGRAVS_ACCUMULATOR
	    Nparticles[i] = nop->u.d.Nparticles[i];
#endif

	    mass[i] = nop->u.d.mass[i];
	    summass += nop->u.d.mass[i];
	    dx[i] = nop->u.d.s[0][i] - pos_x;
	    dy[i] = nop->u.d.s[1][i] - pos_y;
	    dz[i] = nop->u.d.s[2][i] - pos_z;
	    
#ifdef PERIODIC
	    dx[i] = NEAREST(dx[i]);
	    dy[i] = NEAREST(dy[i]);
	    dz[i] = NEAREST(dz[i]);
#endif
	    r2[i] = dx[i] * dx[i] + dy[i] * dy[i] + dz[i] * dz[i];
	    
	    if(r2[i] < r2min)
	      r2min = r2[i];
	    if(r2[i] > r2max)
	      r2max = r2[i];
	    
	  }
	
	  // Flag that we need to be iterating over all gravtypes
	  sG = -1;
	}


      if(no < All.MaxPart)
	{
#ifdef UNEQUALSOFTENINGS
#ifdef ADAPTIVE_GRAVSOFT_FORGAS
	  if(ptype == 0)
	    h = soft;
	  else
	    h = All.ForceSoftening[ptype];

	  if(P[no].Type == 0)
	    {
	      if(h < SphP[no].Hsml)
		h = SphP[no].Hsml;
	    }
	  else
	    {
	      if(h < All.ForceSoftening[P[no].Type])
		h = All.ForceSoftening[P[no].Type];
	    }
#else
	  h = All.ForceSoftening[ptype];
	  if(h < All.ForceSoftening[P[no].Type])
	    h = All.ForceSoftening[P[no].Type];
#endif
#endif
	  no = Nextnode[no];
	}
      else			/* we have an internal node. Need to check opening criterion */
	{
	  if(mode == 1)
	    {
	      if((nop->u.d.bitflags & 3) == 1)	/* if it's a top-level node
						 * which does not contain
						 * local particles we can make
						 * a short-cut 
						 */
		{
		  no = nop->u.d.sibling;
		  continue;
		}
	    }

	  if(All.ErrTolTheta)	/* check Barnes-Hut opening criterion */
	    {
	      if(nop->len * nop->len > r2min * All.ErrTolTheta * All.ErrTolTheta)
		{
		  /* open cell */
		  no = nop->u.d.nextnode;
		  continue;
		}
	    }
	  else			/* check relative opening criterion */
	    {
	      if(summass * nop->len * nop->len > r2min * r2min * aold)
		{
		  /* open cell */
		  no = nop->u.d.nextnode;
		  continue;
		}

	      if(fabs(nop->center[0] - pos_x) < 0.60 * nop->len)
		{
		  if(fabs(nop->center[1] - pos_y) < 0.60 * nop->len)
		    {
		      if(fabs(nop->center[2] - pos_z) < 0.60 * nop->len)
			{
			  no = nop->u.d.nextnode;
			  continue;
			}
		    }
		}
	    }
#ifdef UNEQUALSOFTENINGS
#ifndef ADAPTIVE_GRAVSOFT_FORGAS
	  h = All.ForceSoftening[ptype];
          maxsofttype = (nop->u.d.bitflags >> 2) & 7;
          if(maxsofttype == 7) /* may only occur for zero mass top-level nodes */
            {
              if(mass > 0)
                endrun(988);
              no = nop->u.d.nextnode;
              continue;
            }
          else
            {
              if(h < All.ForceSoftening[maxsofttype])
                {
                  h = All.ForceSoftening[maxsofttype];
                  if(r2max < h * h)
                    {
                      if(((nop->u.d.bitflags >> 5) & 1))	/* bit-5 signals that there are particles of different softening in the node */
                        {
                          no = nop->u.d.nextnode;
                          continue;
                        }
                    }
                }
            }
#else
	  if(ptype == 0)
	    h = soft;
	  else
	    h = All.ForceSoftening[ptype];

	  if(h < nop->maxsoft)
	    {
	      h = nop->maxsoft;
	      if(r2max < h * h)
		{
		  no = nop->u.d.nextnode;
		  continue;
		}
	    }
#endif
#endif

	  no = nop->u.d.sibling;	/* node can be used */

	  if(mode == 1)
	    {
	      if(((nop->u.d.bitflags) & 1))	/* Bit 0 signals that this node belongs to top-level tree */
		continue;
	    }
	}

      pgravtype = TypeToGrav[ptype];

      if(sG > -1) {

	r = sqrt(r2[sG]);
	
	if(r >= h)
	  pot -= (*PotentialFxns[pgravtype][sG])(pmass, mass[sG], h, r, 1);
	else
	  pot += (*PotentialSplines[pgravtype][sG])(pmass, mass[sG], h, r, 1);
	
#ifdef PERIODIC
	pot += mass[sG] * lattice_pot_corr(dx[sG], dy[sG], dz[sG], pgravtype, sG);
#endif
      }
      else {

	// Iterate
	for(i = 0; i < N_GRAVS; ++i) {

	  if(mass[i] == 0.0)
	    continue;

	  r = sqrt(r2[i]);
	
	  if(r >= h) {
#ifdef NGRAVS_ACCUMULATOR
	    pot -= (*PotentialFxns[pgravtype][i])(pmass, mass[i], h, r, Nparticles[i]);
#else
	    pot -= (*PotentialFxns[pgravtype][i])(pmass, mass[i], h, r, 1);
#endif
	  else
	    {
#ifdef NGRAVS_ACCUMULATOR
	      pot += (*PotentialSplines[pgravtype][i])(pmass, mass[i], h, r, Nparticles[i]);
#else
	      pot += (*PotentialSplines[pgravtype][i])(pmass, mass[i], h, r, 1);
#endif
	    }
	  
#ifdef PERIODIC
	  pot += mass[i] * lattice_pot_corr(dx[i], dy[i], dz[i], pgravtype, i);
#endif
	}
      }
    }
	  
  /* store result at the proper place */
  
  if(mode == 0)
    P[target].Potential = pot;
  else
    GravDataResult[target].u.Potential = pot;
#endif
}
	  



#ifdef PMGRID
/*! This function computes the short-range potential when the TreePM
 *  algorithm is used. This potential is the Newtonian potential, modified
 *  by a complementary error function.
 */
void force_treeevaluate_potential_shortrange(int target, int mode)
{

  struct NODE *nop = 0;
  int no, ptype, tabindex;
  double r2[N_GRAVS], dx[N_GRAVS], dy[N_GRAVS], dz[N_GRAVS], mass[N_GRAVS], r, h;
  double pot, pos_x, pos_y, pos_z, aold;
  double eff_dist, fac, rcut, asmth, asmthfac;
  double dxx, dyy, dzz;
#if defined(UNEQUALSOFTENINGS) && !defined(ADAPTIVE_GRAVSOFT_FORGAS)
  int maxsofttype;
#endif
#ifdef ADAPTIVE_GRAVSOFT_FORGAS
  double soft = 0;
#endif

  // ngravs extension
  double summass, r2min, r2max;
  FLOAT pmass;
  int i, pgravtype, sG;
#ifdef NGRAVS_ACCUMULATOR
  long Nparticles[N_GRAVS];
#endif

#ifdef PERIODIC
  double boxsize, boxhalf;
  boxsize = All.BoxSize;
  boxhalf = 0.5 * All.BoxSize;
#endif

  pot = 0;

  if(mode == 0)
    {
      pos_x = P[target].Pos[0];
      pos_y = P[target].Pos[1];
      pos_z = P[target].Pos[2];
      ptype = P[target].Type;
      pmass = P[target].Mass;

      aold = All.ErrTolForceAcc * P[target].OldAcc;
#ifdef ADAPTIVE_GRAVSOFT_FORGAS
      if(ptype == 0)
	soft = SphP[target].Hsml;
#endif
    }
  else
    {
      pos_x = GravDataGet[target].u.Pos[0];
      pos_y = GravDataGet[target].u.Pos[1];
      pos_z = GravDataGet[target].u.Pos[2];
      pmass = GravDataGet[target].Mass;
      ptype = GravDataGet[target].Type;
      aold = All.ErrTolForceAcc * GravDataGet[target].w.OldAcc;
#ifdef ADAPTIVE_GRAVSOFT_FORGAS
      if(ptype == 0)
	soft = GravDataGet[target].Soft;
#endif
    }


  rcut = All.Rcut[0];
  asmth = All.Asmth[0];
#ifdef PLACEHIGHRESREGION
  if(((1 << ptype) & (PLACEHIGHRESREGION)))
    {
      rcut = All.Rcut[1];
      asmth = All.Asmth[1];
    }
#endif
  asmthfac = 0.5 / asmth * (NTAB / 3.0);

#ifndef UNEQUALSOFTENINGS
  h = All.ForceSoftening[ptype];
#endif

  no = All.MaxPart;

  while(no >= 0)
    {
      if(no < All.MaxPart)	/* single particle */
	{
	  /* the index of the node is the index of the particle */
	  /* observe the sign  */
	  sG = TypeToGrav[P[no].Type];
	  mass[sG] = P[no].Mass;	  
	  dx[sG] = P[no].Pos[0] - pos_x;
	  dy[sG] = P[no].Pos[1] - pos_y;
	  dz[sG] = P[no].Pos[2] - pos_z;

#ifdef PERIODIC
	  dx[sG] = NEAREST(dx[sG]);
	  dy[sG] = NEAREST(dy[sG]);
	  dz[sG] = NEAREST(dz[sG]);
#endif
	  r2[sG] = dx[sG] * dx[sG] + dy[sG] * dy[sG] + dz[sG] * dz[sG];
#ifdef NGRAVS_ACCUMULATOR
	  Nparticles[sG] = 1;
#endif
	}
      else
	{
	  if(no >= All.MaxPart + MaxNodes)	/* pseudo particle */
	    {
	      if(mode == 0)
		{
		  Exportflag[DomainTask[no - (All.MaxPart + MaxNodes)]] = 1;
		}
	      no = Nextnode[no - MaxNodes];
	      continue;
	    }

	  r2min = INFINITY;
	  r2max = -INFINITY;
	  summass = 0;

	  nop = &Nodes[no];

	  for(i = 0; i < N_GRAVS; ++i) {

#ifdef NGRAVS_ACCUMULATOR
	    Nparticles[i] = nop->u.d.Nparticles[i];
#endif
	    mass[i] = nop->u.d.mass[i];
	    summass += nop->u.d.mass[i];
	    dx[i] = nop->u.d.s[0][i] - pos_x;
	    dy[i] = nop->u.d.s[1][i] - pos_y;
	    dz[i] = nop->u.d.s[2][i] - pos_z;
	    
#ifdef PERIODIC
	    dx[i] = NEAREST(dx[i]);
	    dy[i] = NEAREST(dy[i]);
	    dz[i] = NEAREST(dz[i]);
#endif
	    r2[i] = dx[i] * dx[i] + dy[i] * dy[i] + dz[i] * dz[i];
	    	    
	    if(r2[i] < r2min)
	      r2min = r2[i];
	    if(r2[i] > r2max)
	      r2max = r2[i];
	  }

	  // Flag to iterate over N_GRAVS
	  sG = -1;
	}

      if(no < All.MaxPart)
	{
#ifdef UNEQUALSOFTENINGS
#ifdef ADAPTIVE_GRAVSOFT_FORGAS
	  if(ptype == 0)
	    h = soft;
	  else
	    h = All.ForceSoftening[ptype];

	  if(P[no].Type == 0)
	    {
	      if(h < SphP[no].Hsml)
		h = SphP[no].Hsml;
	    }
	  else
	    {
	      if(h < All.ForceSoftening[P[no].Type])
		h = All.ForceSoftening[P[no].Type];
	    }
#else
	  h = All.ForceSoftening[ptype];
	  if(h < All.ForceSoftening[P[no].Type])
	    h = All.ForceSoftening[P[no].Type];
#endif
#endif
	  no = Nextnode[no];
	}
      else			/* we have an  internal node. Need to check opening criterion */
	{
	  /* check whether we can stop walking along this branch */
	  if(no >= All.MaxPart + MaxNodes)	/* pseudo particle */
	    {
	      if(mode == 0)
		{
		  Exportflag[DomainTask[no - (All.MaxPart + MaxNodes)]] = 1;
		}
	      no = Nextnode[no - MaxNodes];
	      continue;
	    }

	  if(mode == 1)
	    {
	      if((nop->u.d.bitflags & 3) == 1)	/* if it's a top-level node which does not contain local particles */
		{
		  no = nop->u.d.sibling;
		  continue;
		}
	    }

	  // He has no cutoff check?

	  eff_dist = rcut + 0.5 * nop->len;

	  dxx = nop->center[0] - pos_x;	/* observe the sign ! */
	  dyy = nop->center[1] - pos_y;	/* this vector is -y in my thesis notation */
	  dzz = nop->center[2] - pos_z;
#ifdef PERIODIC
	  dxx = NEAREST(dxx);
	  dyy = NEAREST(dyy);
	  dzz = NEAREST(dzz);
#endif
	  if(dxx < -eff_dist || dxx > eff_dist)
	    {
	      no = nop->u.d.sibling;
	      continue;
	    }

	  if(dyy < -eff_dist || dyy > eff_dist)
	    {
	      no = nop->u.d.sibling;
	      continue;
	    }

	  if(dzz < -eff_dist || dzz > eff_dist)
	    {
	      no = nop->u.d.sibling;
	      continue;
	    }

	  if(All.ErrTolTheta)	/* check Barnes-Hut opening criterion */
	    {
	      if(nop->len * nop->len > r2min * All.ErrTolTheta * All.ErrTolTheta)
		{
		  /* open cell */
		  no = nop->u.d.nextnode;
		  continue;
		}
	    }
	  else			/* check relative opening criterion */
	    {
	      if(summass * nop->len * nop->len > r2min * r2min * aold)
		{
		  /* open cell */
		  no = nop->u.d.nextnode;
		  continue;
		}

	      if(fabs(nop->center[0] - pos_x) < 0.60 * nop->len)
		{
		  if(fabs(nop->center[1] - pos_y) < 0.60 * nop->len)
		    {
		      if(fabs(nop->center[2] - pos_z) < 0.60 * nop->len)
			{
			  no = nop->u.d.nextnode;
			  continue;
			}
		    }
		}
	    }

#ifdef UNEQUALSOFTENINGS
#ifndef ADAPTIVE_GRAVSOFT_FORGAS
	  h = All.ForceSoftening[ptype];
          maxsofttype = (nop->u.d.bitflags >> 2) & 7;
          if(maxsofttype == 7) /* may only occur for zero mass top-level nodes */
            {
              if(summass > 0)
                endrun(989);
              no = nop->u.d.nextnode;
              continue;
            }
          else
            {
              if(h < All.ForceSoftening[maxsofttype])
                {
                  h = All.ForceSoftening[maxsofttype];
                  if(r2max < h * h)
                    {
                      /* bit-5 signals that there are particles of
                       * different softening in the node
                       */
                      if(((nop->u.d.bitflags >> 5) & 1))
                        {
                          no = nop->u.d.nextnode;
                          continue;
                        }
                    }
                }
	    }
#else
	  if(ptype == 0)
	    h = soft;
	  else
	    h = All.ForceSoftening[ptype];

	  if(h < nop->maxsoft)
	    {
	      h = nop->maxsoft;
	      if(r2max < h * h)
		{
		  no = nop->u.d.nextnode;
		  continue;
		}
	    }
#endif
#endif
	  no = nop->u.d.sibling;	/* node can be used */

	  if(mode == 1)
	    {
	      if(((nop->u.d.bitflags) & 1))	/* Bit 0 signals that this node belongs to top-level tree */
		continue;
	    }
	}

      pgravtype = TypeToGrav[ptype];
      
      if(sG > -1) {
	
	r = sqrt(r2[sG]);

	tabindex = (int) (r * asmthfac);

	if(tabindex < NTAB)
	  {
	    fac = shortrange_table_potential[tabindex];

	    if(r >= h)
	      pot -= fac * (*PotentialFxns[pgravtype][sG])(pmass, mass[sG], h, r, 1);
	    else
	      pot += fac * (*PotentialSplines[pgravtype][sG])(pmass, mass[sG], h, r, 1);
	  }
      }
      else {

	for(i = 0; i < N_GRAVS; ++i) {

	  if(mass[i] == 0.0)
	    continue;

	  r = sqrt(r2[i]);

	  tabindex = (int) (r * asmthfac);

	  if(tabindex < NTAB)
	    {
	      fac = shortrange_table_potential[tabindex];

	      if(r >= h) {
#ifdef NGRAVS_ACCUMULATOR
		pot -= fac * (*PotentialFxns[pgravtype][i])(pmass, mass[i], h, r, Nparticles[i]);
#else
		pot -= fac * (*PotentialFxns[pgravtype][i])(pmass, mass[i], h, r, 1);
#endif
	      }
	      else
		{
#ifdef NGRAVS_ACCUMULATOR
		  pot += fac * (*PotentialSplines[pgravtype][i])(pmass, mass[i], h, r, Nparticles[i]);
#else
		  pot += fac * (*PotentialSplines[pgravtype][i])(pmass, mass[i], h, r, 1);
#endif
		}
	    }
	}
	
      }
    }

  /* store result at the proper place */
  if(mode == 0)
    P[target].Potential = pot;
  else
    GravDataResult[target].u.Potential = pot;
}

#endif



/*! This function allocates the memory used for storage of the tree and of
 *  auxiliary arrays needed for tree-walk and link-lists.  Usually,
 *  maxnodes approximately equal to 0.7*maxpart is sufficient to store the
 *  tree for up to maxpart particles.
 */
void force_treeallocate(int maxnodes, int maxpart)
{
  int i;
  size_t bytes;
  double allbytes = 0;
  double u;

  MaxNodes = maxnodes;

  if(!(Nodes_base = malloc(bytes = (MaxNodes + 1) * sizeof(struct NODE))))
    {
      printf("failed to allocate memory for %d tree-nodes (%g MB).\n", MaxNodes, bytes / (1024.0 * 1024.0));
      endrun(3);
    }
  allbytes += bytes;

  if(!(Extnodes_base = malloc(bytes = (MaxNodes + 1) * sizeof(struct extNODE))))
    {
      printf("failed to allocate memory for %d tree-extnodes (%g MB).\n", MaxNodes,
	     bytes / (1024.0 * 1024.0));
      endrun(3);
    }
  allbytes += bytes;

  Nodes = Nodes_base - All.MaxPart;
  Extnodes = Extnodes_base - All.MaxPart;

  if(!(Nextnode = malloc(bytes = (maxpart + MAXTOPNODES) * sizeof(int))))
    {
      printf("Failed to allocate %d spaces for 'Nextnode' array (%g MB)\n", maxpart + MAXTOPNODES,
	     bytes / (1024.0 * 1024.0));
      exit(0);
    }
  allbytes += bytes;

  if(!(Father = malloc(bytes = (maxpart) * sizeof(int))))
    {
      printf("Failed to allocate %d spaces for 'Father' array (%g MB)\n", maxpart, bytes / (1024.0 * 1024.0));
      exit(0);
    }
  allbytes += bytes;

  if(first_flag == 0)
    {
      first_flag = 1;

      if(ThisTask == 0)
	printf("\nAllocated %g MByte for BH-tree. %ld\n\n", allbytes / (1024.0 * 1024.0),
	       sizeof(struct NODE) + sizeof(struct extNODE));

      tabfac = NTAB / 3.0;

      // KC 9/12/15
      // XXX
      // Note concealed Newtonian assumption in the x-space softening tabulations!
      for(i = 0; i < NTAB; i++)
	{
	  u = 3.0 / NTAB * (i + 0.5);
	  shortrange_table[i] = erfc(u) + 2.0 * u / sqrt(M_PI) * exp(-u * u);
	  shortrange_table_potential[i] = erfc(u);
	}
    }
}


/*! This function frees the memory allocated for the tree, i.e. it frees
 *  the space allocated by the function force_treeallocate().
 */
void force_treefree(void)
{
  free(Father);
  free(Nextnode);
  free(Extnodes_base);
  free(Nodes_base);
}




/*! This function does the force computation with direct summation for the
 *  specified particle in the communication buffer. This can be useful for
 *  debugging purposes, in particular for explicit checks of the force
 *  accuracy.
 */
#ifdef FORCETEST
int force_treeevaluate_direct(int target, int mode)
{
  double epsilon;
  double h, h_inv, dx, dy, dz, r, r2, u, fac;
  int i, ptype;
  double pos_x, pos_y, pos_z;
  double acc_x, acc_y, acc_z;
  double pmass;
  
  int pgravtype;
#ifdef PERIODIC
  double fcorr[3];
  //#endif
  //#ifdef PERIODIC
  double boxsize, boxhalf;
 
  boxsize = All.BoxSize;
  boxhalf = 0.5 * All.BoxSize;
#endif

  acc_x = 0;
  acc_y = 0;
  acc_z = 0;

  if(mode == 0)
    {
      pos_x = P[target].Pos[0];
      pos_y = P[target].Pos[1];
      pos_z = P[target].Pos[2];
      ptype = P[target].Type;
      pmass = P[target].Mass;
    }
  else
    {
      pos_x = GravDataGet[target].u.Pos[0];
      pos_y = GravDataGet[target].u.Pos[1];
      pos_z = GravDataGet[target].u.Pos[2];
      pmass = GravDataGet[target].Mass;
      ptype = GravDataGet[target].Type;
    }

  pgravtype = TypeToGrav[ptype];

  for(i = 0; i < NumPart; i++)
    {
      epsilon = dmax(All.ForceSoftening[P[i].Type], All.ForceSoftening[ptype]);

      h = epsilon;
      h_inv = 1 / h;

      dx = P[i].Pos[0] - pos_x;
      dy = P[i].Pos[1] - pos_y;
      dz = P[i].Pos[2] - pos_z;

#ifdef PERIODIC
      while(dx > boxhalf)
	dx -= boxsize;
      while(dy > boxhalf)
	dy -= boxsize;
      while(dz > boxhalf)
	dz -= boxsize;
      while(dx < -boxhalf)
	dx += boxsize;
      while(dy < -boxhalf)
	dy += boxsize;
      while(dz < -boxhalf)
	dz += boxsize;
#endif
      r2 = dx * dx + dy * dy + dz * dz;

      r = sqrt(r2);

      u = r * h_inv;

      if(u >= 1)
	fac = (*AccelFxns[pgravtype][TypeToGrav[P[i].Type]])(pmass, P[i].Mass, h, r, 1);
      else
	fac = (*AccelSplines[pgravtype][TypeToGrav[P[i].Type]])(pmass, P[i].Mass, h, r, 1);
	
      acc_x += dx * fac;
      acc_y += dy * fac;
      acc_z += dz * fac;

#ifdef PERIODIC
      if(u > 1.0e-5)
	{
	  lattice_corr(dx, dy, dz, pgravtype, TypeToGrav[P[i].Type], fcorr);

	  acc_x += P[i].Mass * fcorr[0];
	  acc_y += P[i].Mass * fcorr[1];
	  acc_z += P[i].Mass * fcorr[2];
	}
#endif
    }


  if(mode == 0)
    {
      P[target].GravAccelDirect[0] = acc_x;
      P[target].GravAccelDirect[1] = acc_y;
      P[target].GravAccelDirect[2] = acc_z;
    }
  else
    {
      GravDataResult[target].u.Acc[0] = acc_x;
      GravDataResult[target].u.Acc[1] = acc_y;
      GravDataResult[target].u.Acc[2] = acc_z;
    }


  return NumPart;
}
#endif


/*! This function dumps some of the basic particle data to a file. In case
 *  the tree construction fails, it is called just before the run
 *  terminates with an error message. Examination of the generated file may
 *  then give clues to what caused the problem.
 */
void dump_particles(void)
{
  FILE *fd;
  char buffer[200];
  int i;

  sprintf(buffer, "particles%d.dat", ThisTask);
  fd = fopen(buffer, "w");
  my_fwrite(&NumPart, 1, sizeof(int), fd);

  for(i = 0; i < NumPart; i++)
    my_fwrite(&P[i].Pos[0], 3, sizeof(FLOAT), fd);

  for(i = 0; i < NumPart; i++)
    my_fwrite(&P[i].Vel[0], 3, sizeof(FLOAT), fd);

  for(i = 0; i < NumPart; i++)
    my_fwrite(&P[i].ID, 1, sizeof(int), fd);

  fclose(fd);
}



#ifdef PERIODIC

/*! This function initializes tables with the correction force and the
 *  correction potential due to the periodic images of a point source located
 *  at the origin. 
 *
 *  The corrections specific to the Newtonian interaction 
 *  are obtained by Ewald summation using the original routines given
 *  in Gadget-2 following (Hernquist, Bouchet, Suto, ApJS, 1991, 75, 231) 
 *
 *  Exotic density contributions can be implemented by adjusting the Bertaut procedure
 *  to the particular density at hand: usually just dropping the correction term A in
 *  Eqn. (2.89) (Glasser, Zucker; Theoretical Chemistry:
 *  Advances and Perspectives, 1980).  Care must be taken by the *user* to guarantee
 *  that the routines defined in LatticeForce[][] and LatticeSelf[][] 
 *  give good approximations to the lattice sum of whatever distribution is under
 *  study.  
 *  
 *  The correction fields
 *  are used to obtain the full periodic force if periodic boundaries
 *  combined with the pure tree algorithm are used. For the TreePM
 *  algorithm, these lattice sums are not used.
 *
 *  The correction fields are stored on disk once they are computed. If a
 *  corresponding file is found, they are loaded from disk to speed up the
 *  initialization.  The summation is done in parallel, i.e. the
 *  processors share the work to compute the tables if needed.
 */
void lattice_init(void)
{
  int i, j, k, beg, len, size, n, task, count;
  double x[3], force[3];
  char buf[200];
  FILE *fd;
  
  int l, m;

  if(ThisTask == 0)
    {
      printf("ngravs: initialize Lattice sum correction...\n");
      fflush(stdout);
    }
  
  for(l = 0; l < N_GRAVS; ++l) {
    for(m = 0; m < N_GRAVS; ++m) {
      
      // KC 12/5/14
      // 
#ifdef DOUBLEPRECISION
      sprintf(buf, "lattice_spc_table_%d_dbl_%s.dat", EN, NgravsNames[l][m]);
#else
      sprintf(buf, "lattice_spc_table_%d_%s.dat", EN, NgravsNames[l][m]);
#endif

      if((fd = fopen(buf, "r")))
	{
	  if(ThisTask == 0)
	    {
	      printf("\nngravs: reading Lattice sum tables from file `%s'\n", buf);
	      fflush(stdout);
	    }

	  my_fread(&fcorrx[l][m][0][0][0], sizeof(FLOAT), (EN + 1) * (EN + 1) * (EN + 1), fd);
	  my_fread(&fcorry[l][m][0][0][0], sizeof(FLOAT), (EN + 1) * (EN + 1) * (EN + 1), fd);
	  my_fread(&fcorrz[l][m][0][0][0], sizeof(FLOAT), (EN + 1) * (EN + 1) * (EN + 1), fd);
	  my_fread(&potcorr[l][m][0][0][0], sizeof(FLOAT), (EN + 1) * (EN + 1) * (EN + 1), fd);
	  fclose(fd);
	}
      else
	{
	  if(ThisTask == 0)
	    {
	      printf("\nNo Lattice sum tables in file `%s' found.\nRecomputing them...\n", buf);
	      fflush(stdout);
	    }

	  /* ok, let's recompute things. Actually, we do that in parallel. */

	  size = (EN + 1) * (EN + 1) * (EN + 1) / NTask;


	  beg = ThisTask * size;
	  len = size;
	  if(ThisTask == (NTask - 1))
	    len = (EN + 1) * (EN + 1) * (EN + 1) - beg;

	  for(i = 0, count = 0; i <= EN; i++)
	    for(j = 0; j <= EN; j++)
	      for(k = 0; k <= EN; k++)
		{
		  n = (i * (EN + 1) + j) * (EN + 1) + k;
		  if(n >= beg && n < (beg + len))
		    {
		      if(ThisTask == 0)
			{
			  if((count % (len / 20)) == 0)
			    {
			      printf("%4.1f percent done\n", count / (len / 100.0));
			      fflush(stdout);
			    }
			}

		      x[0] = 0.5 * ((double) i) / EN;
		      x[1] = 0.5 * ((double) j) / EN;
		      x[2] = 0.5 * ((double) k) / EN;

		      (*LatticeForce[l][m])(i, j, k, x, force);

		      fcorrx[l][m][i][j][k] = force[0];
		      fcorry[l][m][i][j][k] = force[1];
		      fcorrz[l][m][i][j][k] = force[2];

		      // KC 12/5/14
		      if(i + j + k == 0)
			potcorr[l][m][i][j][k] = LatticeZero[l][m];
		      else
			potcorr[l][m][i][j][k] = (*LatticePotential[l][m])(x);

		      count++;
		    }
		}

	  for(task = 0; task < NTask; task++)
	    {
	      beg = task * size;
	      len = size;
	      if(task == (NTask - 1))
		len = (EN + 1) * (EN + 1) * (EN + 1) - beg;

#ifdef DOUBLEPRECISION
	      MPI_Bcast(&fcorrx[l][m][0][0][beg], len, MPI_DOUBLE, task, MPI_COMM_WORLD);
	      MPI_Bcast(&fcorry[l][m][0][0][beg], len, MPI_DOUBLE, task, MPI_COMM_WORLD);
	      MPI_Bcast(&fcorrz[l][m][0][0][beg], len, MPI_DOUBLE, task, MPI_COMM_WORLD);
	      MPI_Bcast(&potcorr[l][m][0][0][beg], len, MPI_DOUBLE, task, MPI_COMM_WORLD);
#else
	      MPI_Bcast(&fcorrx[l][m][0][0][beg], len, MPI_FLOAT, task, MPI_COMM_WORLD);
	      MPI_Bcast(&fcorry[l][m][0][0][beg], len, MPI_FLOAT, task, MPI_COMM_WORLD);
	      MPI_Bcast(&fcorrz[l][m][0][0][beg], len, MPI_FLOAT, task, MPI_COMM_WORLD);
	      MPI_Bcast(&potcorr[l][m][0][0][beg], len, MPI_FLOAT, task, MPI_COMM_WORLD);
#endif
	    }

	  if(ThisTask == 0)
	    {
	      printf("\nwriting Lattice sum tables to file `%s'\n", buf);
	      fflush(stdout);

	      if((fd = fopen(buf, "w")))
		{
		  my_fwrite(&fcorrx[l][m][0][0][0], sizeof(FLOAT), (EN + 1) * (EN + 1) * (EN + 1), fd);
		  my_fwrite(&fcorry[l][m][0][0][0], sizeof(FLOAT), (EN + 1) * (EN + 1) * (EN + 1), fd);
		  my_fwrite(&fcorrz[l][m][0][0][0], sizeof(FLOAT), (EN + 1) * (EN + 1) * (EN + 1), fd);
		  my_fwrite(&potcorr[l][m][0][0][0], sizeof(FLOAT), (EN + 1) * (EN + 1) * (EN + 1), fd);
		  fclose(fd);
		}
	    }
	  // KC 12/16/14
	  // Gotta add an MPI_Barrier so that other threads don't run ahead and try to read things
	  MPI_Barrier(MPI_COMM_WORLD);
	}

      fac_intp = 2 * EN / All.BoxSize;

      for(i = 0; i <= EN; i++)
	for(j = 0; j <= EN; j++)
	  for(k = 0; k <= EN; k++)
	    {
	      potcorr[l][m][i][j][k] /= All.BoxSize;
	      fcorrx[l][m][i][j][k] /= All.BoxSize * All.BoxSize;
	      fcorry[l][m][i][j][k] /= All.BoxSize * All.BoxSize;
	      fcorrz[l][m][i][j][k] /= All.BoxSize * All.BoxSize;
	    }

        // Close NGRAVS loops
    }
  }
  if(ThisTask == 0)
    {
      printf("initialization of periodic boundaries finished.\n");
      fflush(stdout);
    }
}


/*! This function looks up the correction force due to the infinite number
 *  of periodic particle/node images. We here use trilinear interpolation
 *  to get it from the precomputed tables, which contain one octant
 *  around the target particle at the origin. The other octants are
 *  obtained from it by exploiting the symmetry properties.
 */
#ifdef FORCETEST
void lattice_corr(double dx, double dy, double dz, int target, int source, double *fper)
{
  int signx, signy, signz;
  int i, j, k;
  double u, v, w;
  double f1, f2, f3, f4, f5, f6, f7, f8;

  if(dx < 0)
    {
      dx = -dx;
      signx = +1;
    }
  else
    signx = -1;

  if(dy < 0)
    {
      dy = -dy;
      signy = +1;
    }
  else
    signy = -1;

  if(dz < 0)
    {
      dz = -dz;
      signz = +1;
    }
  else
    signz = -1;

  u = dx * fac_intp;
  i = (int) u;
  if(i >= EN)
    i = EN - 1;
  u -= i;
  v = dy * fac_intp;
  j = (int) v;
  if(j >= EN)
    j = EN - 1;
  v -= j;
  w = dz * fac_intp;
  k = (int) w;
  if(k >= EN)
    k = EN - 1;
  w -= k;

  f1 = (1 - u) * (1 - v) * (1 - w);
  f2 = (1 - u) * (1 - v) * (w);
  f3 = (1 - u) * (v) * (1 - w);
  f4 = (1 - u) * (v) * (w);
  f5 = (u) * (1 - v) * (1 - w);
  f6 = (u) * (1 - v) * (w);
  f7 = (u) * (v) * (1 - w);
  f8 = (u) * (v) * (w);

  fper[0] = signx * (fcorrx[target][source][i][j][k] * f1 +
		     fcorrx[target][source][i][j][k + 1] * f2 +
		     fcorrx[target][source][i][j + 1][k] * f3 +
		     fcorrx[target][source][i][j + 1][k + 1] * f4 +
		     fcorrx[target][source][i + 1][j][k] * f5 +
		     fcorrx[target][source][i + 1][j][k + 1] * f6 +
		     fcorrx[target][source][i + 1][j + 1][k] * f7 + fcorrx[target][source][i + 1][j + 1][k + 1] * f8);

  fper[1] = signy * (fcorry[target][source][i][j][k] * f1 +
		     fcorry[target][source][i][j][k + 1] * f2 +
		     fcorry[target][source][i][j + 1][k] * f3 +
		     fcorry[target][source][i][j + 1][k + 1] * f4 +
		     fcorry[target][source][i + 1][j][k] * f5 +
		     fcorry[target][source][i + 1][j][k + 1] * f6 +
		     fcorry[target][source][i + 1][j + 1][k] * f7 + fcorry[target][source][i + 1][j + 1][k + 1] * f8);

  fper[2] = signz * (fcorrz[target][source][i][j][k] * f1 +
		     fcorrz[target][source][i][j][k + 1] * f2 +
		     fcorrz[target][source][i][j + 1][k] * f3 +
		     fcorrz[target][source][i][j + 1][k + 1] * f4 +
		     fcorrz[target][source][i + 1][j][k] * f5 +
		     fcorrz[target][source][i + 1][j][k + 1] * f6 +
		     fcorrz[target][source][i + 1][j + 1][k] * f7 + fcorrz[target][source][i + 1][j + 1][k + 1] * f8);
}
#endif


/*! This function looks up the correction potential due to the infinite
 *  number of periodic particle/node images. We here use tri-linear
 *  interpolation to get it from the precomputed table, which contains
 *  one octant around the target particle at the origin. The other
 *  octants are obtained from it by exploiting symmetry properties.
 */
double lattice_pot_corr(double dx, double dy, double dz, int target, int source)
{
  int i, j, k;
  double u, v, w;
  double f1, f2, f3, f4, f5, f6, f7, f8;

  if(dx < 0)
    dx = -dx;

  if(dy < 0)
    dy = -dy;

  if(dz < 0)
    dz = -dz;

  u = dx * fac_intp;
  i = (int) u;
  if(i >= EN)
    i = EN - 1;
  u -= i;
  v = dy * fac_intp;
  j = (int) v;
  if(j >= EN)
    j = EN - 1;
  v -= j;
  w = dz * fac_intp;
  k = (int) w;
  if(k >= EN)
    k = EN - 1;
  w -= k;

  f1 = (1 - u) * (1 - v) * (1 - w);
  f2 = (1 - u) * (1 - v) * (w);
  f3 = (1 - u) * (v) * (1 - w);
  f4 = (1 - u) * (v) * (w);
  f5 = (u) * (1 - v) * (1 - w);
  f6 = (u) * (1 - v) * (w);
  f7 = (u) * (v) * (1 - w);
  f8 = (u) * (v) * (w);

  return potcorr[target][source][i][j][k] * f1 +
    potcorr[target][source][i][j][k + 1] * f2 +
    potcorr[target][source][i][j + 1][k] * f3 +
    potcorr[target][source][i][j + 1][k + 1] * f4 +
    potcorr[target][source][i + 1][j][k] * f5 +
    potcorr[target][source][i + 1][j][k + 1] * f6 + potcorr[target][source][i + 1][j + 1][k] * f7 + potcorr[target][source][i + 1][j + 1][k + 1] * f8;
}



/*! This function computes the potential correction term by means of Ewald
 *  summation.
 */
double ewald_psi(double x[3])
{
  double alpha, psi;
  double r, sum1, sum2, hdotx;
  double dx[3];
  int i, n[3], h[3], h2;

  alpha = 2.0;

  for(n[0] = -4, sum1 = 0; n[0] <= 4; n[0]++)
    for(n[1] = -4; n[1] <= 4; n[1]++)
      for(n[2] = -4; n[2] <= 4; n[2]++)
	{
	  for(i = 0; i < 3; i++)
	    dx[i] = x[i] - n[i];

	  r = sqrt(dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2]);
	  sum1 += erfc(alpha * r) / r;
	}

  for(h[0] = -4, sum2 = 0; h[0] <= 4; h[0]++)
    for(h[1] = -4; h[1] <= 4; h[1]++)
      for(h[2] = -4; h[2] <= 4; h[2]++)
	{
	  hdotx = x[0] * h[0] + x[1] * h[1] + x[2] * h[2];
	  h2 = h[0] * h[0] + h[1] * h[1] + h[2] * h[2];
	  if(h2 > 0)
	    sum2 += 1 / (M_PI * h2) * exp(-M_PI * M_PI * h2 / (alpha * alpha)) * cos(2 * M_PI * hdotx);
	}

  r = sqrt(x[0] * x[0] + x[1] * x[1] + x[2] * x[2]);

  psi = M_PI / (alpha * alpha) - sum1 - sum2 + 1 / r;

  return psi;
}


/*! This function computes the force correction term (difference between full
 *  force of infinite lattice and nearest image) by Ewald summation.
 */
void ewald_force(int iii, int jjj, int kkk, double x[3], double force[3])
{
  double alpha, r2;
  double r, val, hdotx, dx[3];
  int i, h[3], n[3], h2;

  alpha = 2.0;

  for(i = 0; i < 3; i++)
    force[i] = 0;

  if(iii == 0 && jjj == 0 && kkk == 0)
    return;

  r2 = x[0] * x[0] + x[1] * x[1] + x[2] * x[2];

  for(i = 0; i < 3; i++)
    force[i] += x[i] / (r2 * sqrt(r2));

  // KC 12/4/14
  // Looks like this takes the first four images out in position space in each direction (so 
  // bracketing by 8 overall)
  for(n[0] = -4; n[0] <= 4; n[0]++)
    for(n[1] = -4; n[1] <= 4; n[1]++)
      for(n[2] = -4; n[2] <= 4; n[2]++)
	{
	  for(i = 0; i < 3; i++)
	    dx[i] = x[i] - n[i];

	  r = sqrt(dx[0] * dx[0] + dx[1] * dx[1] + dx[2] * dx[2]);

	  val = erfc(alpha * r) + 2 * alpha * r / sqrt(M_PI) * exp(-alpha * alpha * r * r);

	  for(i = 0; i < 3; i++)
	    force[i] -= dx[i] / (r * r * r) * val;
	}

  // KC 12/4/14
  // Looks like this takes the first four images in momentum space in each diretion (again
  // bracketing by 8 overall)
  for(h[0] = -4; h[0] <= 4; h[0]++)
    for(h[1] = -4; h[1] <= 4; h[1]++)
      for(h[2] = -4; h[2] <= 4; h[2]++)
	{
	  hdotx = x[0] * h[0] + x[1] * h[1] + x[2] * h[2];
	  h2 = h[0] * h[0] + h[1] * h[1] + h[2] * h[2];

	  if(h2 > 0)
	    {
	      val = 2.0 / ((double) h2) * exp(-M_PI * M_PI * h2 / (alpha * alpha)) * sin(2 * M_PI * hdotx);

	      for(i = 0; i < 3; i++)
		force[i] -= h[i] * val;
	    }
	}
}

#endif
