/*
Martin Fuller, March 2015

This program counts minimal triangle graphs (MTGs) [1] for the sequences https://oeis.org/A000080 and https://oeis.org/A048781 [2].
The method is similar to that described by R. Bowen in [1] in that graphs are grown from smaller MTGs and then checked for minimality and isomorphism.
Run-time was dominated by isomorphism checking; therefore the following changes were made to reduce the number of isomorphic graphs.

0. Triangles are added only at nodes that are nonisomorphic in the source graph.
   The program uses node orbits when it adds triangles at a single node.
   The program uses the generators of the automorphism group when it adds triangles across 2 nodes.

1. Triangles are added in 3 stages so that each output graph is constructed in fewer ways and there are fewer isomorphs to reject.
   The first stage creates all MTGs up to n nodes that have no "hanging" triangles, that is no triangles have 1 node connected
   to the rest of the graph and 2 nodes that are in no other triangle.  This stage starts from and includes the MTG with one triangle.
   There are 4 addition methods which suffice to make all such MTGs :

   i. A chain of 2+ triangles connecting a node to itself.
   ii. A stack of 1+ hanging triangles followed by a chain of 2+ triangles connecting one of the last 2 nodes to itself.
   iii. A stack of 1+ hanging triangles followed by a chain of 1+ triangles connecting the last 2 nodes.
   iv. A chain of 1+ triangles connecting 2 distinct nodes.

   Example chain of 3 triangles: o = old node in the source graph; + = new node to be added
      o---+---+---o
       \ / \ / \ /
        +   +   +
   A stack is similar except only one end is attached the the source graph and the other end is a new node.

   The initial single triangle only accepts addition iv, as otherwise the initial triangle itself would be left hanging.
   Afterwards all 4 additions are used.

2. The second stage adds hanging triangles to nodes that are not uniquely identified at the end of stage 1; that is, their orbit contains 2+ nodes.
   (If the program is required to output all graphs rather than count them, then this stage is extended to include all nodes and stage 3 is ignored.)
   Each hanging triangle adds 2 new nodes, therefore this stage is only used if the number of nodes in the source graph has the required parity.
   The program repeatedly adds a stack of 1+ triangles to a non-hanging node.  The last 2 new nodes are marked as hanging nodes.
   Branching structures are generated when the program adds a stack to the middle of a previous stack.
   The program partially orders nodes by their orbit in the graph at the end of stage 1, or by the stage 1 node to which they are joined.  It then
   adds stacks in non-decreasing order of these orbits.
   If the stage 1 graph is a single triangle then 2 of its nodes are marked as hanging nodes so the source triangle is forced to be at one end of
   the generated graph.

3. The third stage counts the number of ways to add hanging triangles to nodes that are uniquely identified at the end of stage 1.
   The graphs are not generated.  For a single node the count is sequence https://oeis.org/A003080.  Counts for 2+ nodes are calculated by convolution.

Nauty [3] is used for graph canonisation, node orbits and generators of the automorphism group.

Isomorphs are rejected using a "canonical construction path" [4]  The reduction process is the reverse of the additions described in stages 1-3.
  A hanging triangle is removed first if there are any.
  Otherwise any triangle is removed that leaves a smaller MTG.
In both cases, this is followed by repeatedly removing all triangles that become hanging as a result.
The program searches for the first canonically-labelled node of degree 2 which is part of the above type of triangle.
The current construction path is rejected unless the most recent addition has a node in the same orbit as the node selected above.
In this way each isomorphism class is accepted only once without needing to store the graphs.

The program was compiled using the following command line:
   gcc -W -Wall -std=c99 -O4 -Inauty25r9 a000080.c nauty25r9/nauty1.a -o a000080

Total run-time up to n=25 was approximately 24 hours.  Memory use was a few megabytes.
The rate of production fell as n increased: 54800 MTG/sec at n=15, 33600 MTG/sec at n=20, 23200 MTG/sec at n=25.
Stage 1 required most of the run-time despite generating only a fraction of the MTGs (~1/5 at n=20).

References:
[1] R. Bowen, The generation of minimal triangle graphs, Math. Comp. 21 (1967), 248-250
[2] OEIS Foundation Inc. (2011), The On-Line Encyclopedia of Integer Sequences, http://oeis.org.
[3] McKay, B.D. and Piperno, A., Practical Graph Isomorphism, II, J. Symbolic Computation, 60 (2014), pp. 94-112, http://dx.doi.org/10.1016/j.jsc.2013.09.003.
[4] B.D.McKay, Isomorph-free exhaustive generation, J. Algorithms, 26 (1998) 306-324.

$ for ((n=3; n<=25; n++)) do ( ./a000080 $n ) done
>A a000080 3
3 1 1
>Z a000080(3)=1 a048781(3)=1 in 0.00 sec
>A a000080 4
4 1 2
>Z a000080(4)=1 a048781(4)=2 in 0.00 sec
>A a000080 5
5 2 5
>Z a000080(5)=2 a048781(5)=5 in 0.00 sec
>A a000080 6
6 4 13
>Z a000080(6)=4 a048781(6)=13 in 0.00 sec
>A a000080 7
7 9 35
>Z a000080(7)=9 a048781(7)=35 in 0.00 sec
>A a000080 8
8 19 86
>Z a000080(8)=19 a048781(8)=86 in 0.00 sec
>A a000080 9
9 48 250
>Z a000080(9)=48 a048781(9)=250 in 0.00 sec
>A a000080 10
10 117 687
>Z a000080(10)=117 a048781(10)=687 in 0.00 sec
>A a000080 11
11 307 2001
>Z a000080(11)=307 a048781(11)=2001 in 0.00 sec
>A a000080 12
12 821 5897
>Z a000080(12)=821 a048781(12)=5897 in 0.01 sec
>A a000080 13
13 2277 17847
>Z a000080(13)=2277 a048781(13)=17847 in 0.03 sec
>A a000080 14
14 6437 54683
>Z a000080(14)=6437 a048781(14)=54683 in 0.11 sec
>A a000080 15
15 18634 170578
>Z a000080(15)=18634 a048781(15)=170578 in 0.34 sec
>A a000080 16
16 54775 537468
>Z a000080(16)=54775 a048781(16)=537468 in 1.14 sec
>A a000080 17
17 163703 1714312
>Z a000080(17)=163703 a048781(17)=1714312 in 3.68 sec
>A a000080 18
18 495529 5516392
>Z a000080(18)=495529 a048781(18)=5516392 in 12.15 sec
>A a000080 19
19 1518706 17910741
>Z a000080(19)=1518706 a048781(19)=17910741 in 40.51 sec
>A a000080 20
20 4703848 58588966
>Z a000080(20)=4703848 a048781(20)=58588966 in 136.17 sec
>A a000080 21
21 14714754 193035920
>Z a000080(21)=14714754 a048781(21)=193035920 in 476.16 sec
>A a000080 22
22 46444979 640132185
>Z a000080(22)=46444979 a048781(22)=640132185 in 1657.24 sec
>A a000080 23
23 147832051 2135841554
>Z a000080(23)=147832051 a048781(23)=2135841554 in 5857.21 sec
>A a000080 24
24 474229588 7167550473
>Z a000080(24)=474229588 a048781(24)=7167550473 in 19353.90 sec
>A a000080 25
25 1532565644 24186506740
>Z a000080(25)=1532565644 a048781(25)=24186506740 in 65898.77 sec
*/

#define USAGE \
"a000080 [-tgq] n"

#define HELPTEXT \
" Calculate A000080(n)\n\
\n\
    n   the number of vertices\n\
\n\
    -t  Output triples\n\
    -g  Output graphs\n\
    -q  Suppress auxiliary information\n"

#define MAXN WORDSIZE
#define ONE_WORD_SETS

#include "gtools.h"
#include "naututil.h" /* For setsize() */
#include <limits.h>

typedef unsigned char nsize;
#define NSIZE_FORMAT "%hhu"

typedef struct
{
  nsize n;
  boolean quiet;
  boolean outputCount;
  boolean outputTriples;
  boolean outputGraphs;
  boolean generateAllGraphs;
} argsdata;

TLS_ATTR unsigned long a000080;
TLS_ATTR unsigned long a048781;
TLS_ATTR nsize numGeneratorsGlobal;
TLS_ATTR nsize generatorsGlobal[MAXN * MAXN];
DEFAULTOPTIONS_GRAPH(optionsCanon);
DEFAULTOPTIONS_GRAPH(optionsCanonPlusGenerators);
const unsigned long a003080[] = { 1, 1, 2, 5, 13, 37, 111, 345, 1105, 3624, 12099, 41000, 140647, 487440, 1704115, 6002600, 21282235, 75890812, 272000538, 979310627,
    3540297130u, 12845634348u, 46764904745u, 170767429511u, 625314778963u, 2295635155206u, 8447553316546u };
#define A003080_LEN ((unsigned int)(sizeof(a003080) / sizeof(a003080[0])))
#define CACTI_FORESTS_NODES_LEN 32
unsigned long cactiForests[CACTI_FORESTS_NODES_LEN][A003080_LEN]; /* number of triangular cacti forests (A003080) on ordered root nodes [#ordered root nodes][#triangles] */
unsigned int cactiForestsTrianglesLen[CACTI_FORESTS_NODES_LEN]; /* #triangles that causes integer overflow in above array [#ordered root nodes] */

void initCactiForests()
{
  /* No nodes, no triangles => 1 empty forest. No nodes, some triangles => no forests */
  cactiForests[0][0] = 1;
  for (unsigned int j = 1; j < A003080_LEN; j++)
    cactiForests[0][j] = 0;
  cactiForestsTrianglesLen[0] = A003080_LEN;

  /* One node => A003080 */
  for (unsigned int j = 0; j < A003080_LEN; j++)
    cactiForests[1][j] = a003080[j];
  cactiForestsTrianglesLen[1] = A003080_LEN;

  /* 2+ ordered nodes => convolution of A003080 for first node and previous row for remaining nodes */
  unsigned int len = A003080_LEN;
  for (unsigned int i = 2; i < CACTI_FORESTS_NODES_LEN; i++)
  {
    for (unsigned int j = 0; j < len; j++)
    {
      unsigned long count = 0;
      for (unsigned int k = 0; k <= j; k++)
      {
        unsigned long mul = a003080[k] * cactiForests[i-1][j-k];
        if ((count > ULONG_MAX - mul) || (mul / a003080[k] != cactiForests[i-1][j-k]))
          len = j;
        count += mul;
      }
      cactiForests[i][j] = count;
    }
    cactiForestsTrianglesLen[i] = len;
  }
}

void storeGenerators(int count, int *perm, int *orbits, int numorbits, int stabvertex, int n)
{
  (void)count; (void)orbits; (void)numorbits; (void)stabvertex;

  nsize *g = generatorsGlobal + (numGeneratorsGlobal++) * MAXN;
  int *p = perm;
  for (int i = 0; i < n; i++)
    *g++ = *p++;
}

void initOptions()
{
  optionsCanon.getcanon = TRUE;
  optionsCanonPlusGenerators.getcanon = TRUE;
  optionsCanonPlusGenerators.userautomproc = storeGenerators;
}

argsdata
parseArgs(int argc, char *argv[])
{
  argsdata args;
  int argnum = 0;
  boolean badargs = FALSE;

  args.quiet = FALSE;
  args.outputTriples = FALSE;
  args.outputGraphs = FALSE;

  for (int j = 1; !badargs && j < argc; ++j)
  {
    char *arg = argv[j];
    if (arg[0] == '-' && arg[1] != '\0')
    {
      ++arg;
      while (*arg != '\0')
      {
        char sw = *arg++;
             SWBOOLEAN('t', args.outputTriples)
        else SWBOOLEAN('g', args.outputGraphs)
        else SWBOOLEAN('q', args.quiet)
        else badargs = TRUE;
      }
    }
    else
    {
      ++argnum;
      if (argnum == 1) { args.n = (nsize)atoi(arg); badargs = (args.n == 0); }
      else badargs = TRUE;
    }
  }

  args.outputCount = !(args.outputTriples || args.outputGraphs);
  args.generateAllGraphs = (args.outputTriples || args.outputGraphs);

  if (badargs || argnum == 0)
  {
      fprintf(stderr, ">E Usage: %s\n", USAGE);
      GETHELP;
      exit(1);
  }

  if ((args.n < 3) || (args.n > MAXN))
  {
    fprintf(stderr, ">E n=%d must be between 3 and %d inclusive\n", args.n, MAXN);
    exit(1);
  }

  if (!args.quiet)
  {
    fprintf(stderr, ">A a000080");
    if (args.outputTriples || args.outputGraphs)
    {
      fprintf(stderr, " -");
      if (args.outputTriples) fprintf(stderr, "t");
      if (args.outputGraphs) fprintf(stderr, "g");
    }
    fprintf(stderr, " %d\n", args.n);
  }

  return args;
}

void writeTriples(FILE *file, nsize numTriples, const nsize *triples, const char *suffix)
{
  const nsize *t = triples;
  for (nsize i = 0; i < numTriples; i++, t += 3)
    fprintf(file, "(%d,%d,%d)", t[0], t[1], t[2]);
  fprintf(file, suffix);
}

#if  (MAXM==1) && defined(ONE_WORD_SETS)
#define INTERSECTS(s1,s2,m) ((*(s1) & *(s2)) != 0)
#else
#error Not yet implemented
#endif
#define ISARC(g,v,w,m) ISELEMENT(GRAPHROW(g,v,m),w)

boolean isConnectedExcluding(nsize numTriples, const nsize *triples, nsize node1, nsize node2, nsize excludedTriple)
{
  set s1[MAXM], s2[MAXM], c1[MAXM], c2[MAXM];
  EMPTYSET(s1, MAXM);
  EMPTYSET(s2, MAXM);
  ADDELEMENT(s1, node1);
  ADDELEMENT(s2, node2);
  while (TRUE)
  {
    memcpy(c1, s1, sizeof(s1));
    memcpy(c2, s2, sizeof(s2));
    const nsize *t = triples;
    for (int i = 0; i < numTriples; i++, t += 3)
      if (i != excludedTriple)
      {
        if (ISELEMENT(s1, t[0]) || ISELEMENT(s1, t[1]) || ISELEMENT(s1, t[2]))
        {
          ADDELEMENT(s1, t[0]);
          ADDELEMENT(s1, t[1]);
          ADDELEMENT(s1, t[2]);
        }
        if (ISELEMENT(s2, t[0]) || ISELEMENT(s2, t[1]) || ISELEMENT(s2, t[2]))
        {
          ADDELEMENT(s2, t[0]);
          ADDELEMENT(s2, t[1]);
          ADDELEMENT(s2, t[2]);
        }
      }
    if (INTERSECTS(s1,s2,MAXM))
      return TRUE;
    if ((memcmp(c1, s1, sizeof(c1)) == 0) || (memcmp(c2, s2, sizeof(c2)) == 0))
      return FALSE;
  }
}

boolean isMinimal(nsize numTriples, const nsize *triples)
{
  const nsize *t = triples;
  for (nsize i = 0; i < numTriples; i++, t += 3)
    if (isConnectedExcluding(numTriples, triples, t[0], t[2], i)
      && (((i > 0) && (t[1] <= t[-1])) || isConnectedExcluding(numTriples, triples, t[0], t[1], i)))
      return FALSE;
  return TRUE;
}

const nsize* find(const nsize *array, nsize len, nsize value)
{
  for (nsize i = 0; i < len; i++, array++)
    if (*array == value)
      return array;
  return NULL;
}

boolean canRemoveTriple2(nsize numTriples, const nsize *triples, nsize index)
{
  const nsize *t = triples + index * 3;
  return !find(t + 3, (numTriples - index - 1) * 3, t[2])
      || isConnectedExcluding(numTriples, triples, t[0], t[2], index);
}

boolean canRemoveTriple1(nsize numTriples, const nsize *triples, nsize index)
{
  const nsize *t = triples + index * 3;
  return ((index > 0) && (t[1] <= t[-1]))
      || !find(t + 3, (numTriples - index - 1) * 3, t[1])
      || isConnectedExcluding(numTriples, triples, t[0], t[1], index);
}

boolean canRemoveNode(nsize numTriples, const nsize *triples, nsize nodeIndex)
{
  const nsize *t = find(triples, numTriples * 3, nodeIndex);
  if (find(t + 1, (triples + numTriples * 3 - 1) - t, nodeIndex))
    return FALSE;
  nsize index = (t - triples) / 3;
  switch ((t - triples) % 3)
  {
    case 2:
      return canRemoveTriple1(numTriples, triples, index);
    case 1:
      return canRemoveTriple2(numTriples, triples, index);
    default:
      return canRemoveTriple1(numTriples, triples, index) && canRemoveTriple2(numTriples, triples, index);
  }
}

/*
Find the first hanging node in the canonical labelling.
Return true if one of the recently added triples is in the same orbit.
*/
boolean isCanonicalLastHanging(nsize numTriples, nsize prevNumTriples, const nsize *triples, nsize n, const int *lab, const int *orbits, const set *hangingNodes, const int *basedOrbits)
{
  for (int k = 0; k < n; k++)
  {
    int i = lab[k];
    if (ISELEMENT(hangingNodes, i) && (basedOrbits[i] == basedOrbits[n-1]))
    {
      const nsize *t = triples + prevNumTriples * 3;
      for (int j = prevNumTriples; j < numTriples; j++, t += 3)
        if ((orbits[i] == orbits[t[0]]) || (orbits[i] == orbits[t[1]]) || (orbits[i] == orbits[t[2]]))
          return TRUE;
      return FALSE;
    }
  }
  fprintf(stderr, ">E No hanging triples in ");
  writeTriples(stderr, numTriples, triples, "\n");
  return FALSE;
}

/*
Find the first node in the canonical labelling that can be removed by removing 1 triple while leaving a MTS.
Return true if one of the recently added triples is in the same orbit.
*/
boolean isCanonicalLastBased(nsize numTriples, nsize prevNumTriples, const nsize *triples, nsize n, const int *lab, const int *orbits)
{
  set tested[MAXM];
  EMPTYSET(tested, MAXM);
  for (int k = 0; k < n; k++)
  {
    int i = lab[k];
    if (!ISELEMENT(tested, orbits[i]))
    {
      if (canRemoveNode(numTriples, triples, orbits[i]))
      {
        const nsize *t = triples + prevNumTriples * 3;
        for (nsize j = prevNumTriples; j < numTriples; j++, t += 3)
          if ((orbits[i] == orbits[t[0]]) || (orbits[i] == orbits[t[1]]) || (orbits[i] == orbits[t[2]]))
            return TRUE;
        return FALSE;
      }
      ADDELEMENT(tested, orbits[i]);
    }
  }
  fprintf(stderr, ">E No removable nodes in ");
  writeTriples(stderr, numTriples, triples, "\n");
  return FALSE;
}

void setTriple(nsize *newTriple, nsize a, nsize b, nsize c)
{
  *newTriple++ = a;
  *newTriple++ = b;
  *newTriple++ = c;
}

void setTripleGraph(nsize n, graph *g, nsize prevN, const graph *prevG, nsize numNewTriples, const nsize *newTriples)
{
  memcpy(g, prevG, MAXM * (size_t)prevN * sizeof(g[0]));
  memset(g + MAXM * (size_t)prevN, 0, MAXM * (size_t)(n - prevN) * sizeof(g[0]));
  const nsize *t = newTriples;
  for (int i = 0; i < numNewTriples; i++, t += 3)
  {
    ADDONEEDGE(g, t[0], t[1], MAXM);
    ADDONEEDGE(g, t[0], t[2], MAXM);
    ADDONEEDGE(g, t[1], t[2], MAXM);
  }
}

void addGeneratedEdges(graph *g, int v, int w, nsize numGenerators, const nsize *generators)
{
  if (!ISARC(g, v, w, MAXM))
  {
    ADDONEEDGE(g, v, w, MAXM);
    const nsize *gen = generators;
    for (nsize i = 0; i < numGenerators; i++, gen += MAXN)
      addGeneratedEdges(g, gen[v], gen[w], numGenerators, generators);
  }
}

void output(nsize numTriples, const nsize *triples, const graph *g, const argsdata *args)
{
  if (args->outputTriples)
    writeTriples(stdout, numTriples, triples, "\n");
  if (args->outputGraphs)
    writeg6(stdout, (graph*)g, MAXM, args->n);
}

/* Number of triangular cactus forests with ordered root nodes.  With 1 node it is A003080. More nodes were calculated by convolution */
unsigned long countWithCacti(unsigned int numTriangles, unsigned int numNodes)
{
  if ((numNodes >= CACTI_FORESTS_NODES_LEN) || (numTriangles >= cactiForestsTrianglesLen[numNodes]))
  {
    fprintf(stderr, ">E countWithCacti(%d, %d) too large\n", numTriangles, numNodes);
    exit(1);
  }
  return cactiForests[numNodes][numTriangles];
}

void addHangingTriples(nsize numTriples, nsize triples[3 * MAXN], const graph *g, const int *orbits, const set hangingNodes[MAXM], int basedOrbits[MAXN], int minBasedOrbit, const set *isUniqueNode, const argsdata *args);

void addStack(nsize numTriples, nsize triples[3 * MAXN], nsize prevNumTriples, nsize prevN, const graph *prevG, const set prevHangingNodes[MAXM], int basedOrbits[MAXN], const set *isUniqueNode, const argsdata *args)
{
  nsize n = triples[3 * numTriples - 1] + 1;
  graph g[MAXM * n];
  int orbits[n];
  set hangingNodes[MAXM];
  int newBasedOrbit = basedOrbits[triples[prevNumTriples * 3]];
  {
    int lab[n], ptn[n];
    graph canon[MAXM * n];
    statsblk stats;
    setTripleGraph(n, g, prevN, prevG, numTriples - prevNumTriples, triples + prevNumTriples * 3);
    memcpy(hangingNodes, prevHangingNodes, MAXM * sizeof(hangingNodes[0]));
    ADDELEMENT(hangingNodes, triples[numTriples * 3 - 2]);
    ADDELEMENT(hangingNodes, triples[numTriples * 3 - 1]);
    for (nsize i = prevN; i < n; i++)
      basedOrbits[i] = newBasedOrbit;
    densenauty(g, lab, ptn, orbits, &optionsCanon, &stats, MAXM, n, canon);
    if (!isCanonicalLastHanging(numTriples, prevNumTriples, triples, n, lab, orbits, hangingNodes, basedOrbits))
      return;
  }
  addHangingTriples(numTriples, triples, g, orbits, hangingNodes, basedOrbits, newBasedOrbit, isUniqueNode, args);
}

void addHangingTriples(nsize numTriples, nsize triples[3 * MAXN], const graph *g, const int *orbits, const set hangingNodes[MAXM], int basedOrbits[MAXN], int minBasedOrbit, const set *isUniqueNode, const argsdata *args)
{
  nsize n = triples[3 * numTriples - 1] + 1;

  if (n == args->n)
  {
    output(numTriples, triples, g, args);
    a000080++;
    a048781 += numTriples;
    return;
  }

  if (!args->generateAllGraphs)
  {
    /* Count cacti on nodes that were in an orbit by themselves after the based triples */
    unsigned long count = countWithCacti((args->n - n) / 2, setsize((set*)isUniqueNode, MAXM));
    a000080 += count;
    a048781 += count * (numTriples + (args->n - n) / 2);
  }

  /* Add a stack of hanging triples (i, n, n+1) (n, n+2, n+3) etc once in each current orbit on nodes that were in a 2+ orbit after the based triples */
  for (int k = 1, max = (args->n - n) / 2; k <= max; k++)
  {
    setTriple(triples + (numTriples + k - 1) * 3, n + (k - 2) * 2, n + (k - 1) * 2, n + (k - 1) * 2 + 1);
    for (int i = 0; i < n; i++)
      if ((orbits[i] == i) && !ISELEMENT(hangingNodes, i) && (basedOrbits[i] >= minBasedOrbit) && (args->generateAllGraphs || !ISELEMENT(isUniqueNode, i)))
      {
        triples[numTriples * 3] = i;
        addStack(numTriples + k, triples, numTriples, n, g, hangingNodes, basedOrbits, isUniqueNode, args);
      }
  }
}

void hangingTriplesSecond(nsize numTriples, nsize triples[3 * MAXN], const graph *g, const int *orbits, const argsdata *args)
{
  nsize n = triples[3 * numTriples - 1] + 1;
  set hangingNodes[MAXM];
  int basedOrbits[MAXN];
  set isUniqueNode[MAXM];

  EMPTYSET(hangingNodes, MAXM);
  if (numTriples == 1)
  {
    ADDELEMENT(hangingNodes, 1);
    ADDELEMENT(hangingNodes, 2);
  }
  memcpy(basedOrbits, orbits, n * sizeof(orbits[0]));
  EMPTYSET(isUniqueNode, MAXM);
  for (int i = 0; i < n; i++)
    if (orbits[i] == i)
      ADDELEMENT(isUniqueNode, i);
    else
      DELELEMENT(isUniqueNode, orbits[i]);

  addHangingTriples(numTriples, triples, g, orbits, hangingNodes, basedOrbits, 0, isUniqueNode, args);
}

void basedTriplesFirst(nsize numTriples, nsize triples[3 * MAXN], nsize prevNumTriples, nsize prevN, const graph *prevG, const argsdata *args);

/* Add a chain of triples from v to w of each length: (v, n, n+1), (n, n+2, n+3) ... (w, n+2k-3, n+2k-1) */
void addChain(nsize numTriples, nsize triples[3 * MAXN], nsize v, nsize w, nsize prevNumTriples, nsize prevN, const graph *prevG, const argsdata *args)
{
  nsize n = triples[3 * numTriples - 1] + 1;
  if (v != w)
  {
    setTriple(triples + numTriples * 3, v, w, n);
    basedTriplesFirst(numTriples + 1, triples, prevNumTriples, prevN, prevG, args);
  }
  for (int k = 2, max = (args->n - n + 1) / 2; k <= max; k++)
  {
    setTriple(triples + (numTriples + k - 2) * 3, (k == 2) ? v : (n + (k - 3) * 2), n + (k - 2) * 2, n + (k - 2) * 2 + 1);
    setTriple(triples + (numTriples + k - 1) * 3, w, n + (k - 2) * 2, n + (k - 1) * 2);
    basedTriplesFirst(numTriples + k, triples, prevNumTriples, prevN, prevG, args);
  }
}

void basedTriplesFirst(nsize numTriples, nsize triples[3 * MAXN], nsize prevNumTriples, nsize prevN, const graph *prevG, const argsdata *args)
{
  nsize n = triples[3 * numTriples - 1] + 1;
  graph g[MAXM * n], temp[MAXM * n];
  int lab[n], ptn[n], orbits[n];
  statsblk stats;

  setTripleGraph(n, g, prevN, prevG, numTriples - prevNumTriples, triples + prevNumTriples * 3);
  numGeneratorsGlobal = 0;
  densenauty(g, lab, ptn, orbits, &optionsCanonPlusGenerators, &stats, MAXM, n, temp);
  if (!isCanonicalLastBased(numTriples, prevNumTriples, triples, n, lab, orbits))
    return;

  if (n == args->n)
  {
    output(numTriples, triples, g, args);
    a000080++;
    a048781 += numTriples;
    return;
  }

  nsize numGenerators = numGeneratorsGlobal, generators[numGenerators * MAXN];
  memcpy(generators, generatorsGlobal, sizeof(generators[0]) * numGenerators * MAXN);

  if ((args->n - n) % 2 == 0)
    hangingTriplesSecond(numTriples, triples, g, orbits, args);

  if (numTriples > 1)
  {
    /* i. A chain of 2+ triangles connecting a node to itself. */
    /* ii. A stack of 1+ hanging triangles followed by a chain of 2+ triangles connecting one of the last 2 nodes to itself. */
    /* iii. A stack of 1+ hanging triangles followed by a chain of 1+ triangles connecting the last 2 nodes. */
    for (int i = 0; i < n; i++)
      if (orbits[i] == i)
      {
        addChain(numTriples, triples, i, i, numTriples, n, g, args);
        for (int k = 1, max = (args->n - n - 1) / 2; k <= max; k++)
        {
          nsize v = n + (k - 1) * 2;
          setTriple(triples + (numTriples + k - 1) * 3, (k == 1) ? i : (v - 2), v, v + 1);
          addChain(numTriples + k, triples, v, v, numTriples, n, g, args);
          addChain(numTriples + k, triples, v, v + 1, numTriples, n, g, args);
        }
      }
  }

  /* iv. A chain of 1+ triangles connecting 2 distinct nodes. */
  EMPTYGRAPH(temp, MAXM, n);
  for (int i = 0; i < n; i++)
    for (int j = i + 1; j < n; j++)
      if (!ISARC(temp, i, j, m))
      {
        setTriple(triples + numTriples * 3, i, j, n);
        if (isMinimal(numTriples + 1, triples))
          addChain(numTriples, triples, i, j, numTriples, n, g, args);
        addGeneratedEdges(temp, i, j, numGenerators, generators);
      }
}

void calculate(const argsdata *args)
{
  nsize triples[3 * MAXN];
  setTriple(triples, 0, 1, 2);
  a000080 = a048781 = 0;
  basedTriplesFirst(1, triples, 0, 0, NULL, args);
}

int
main(int argc, char *argv[])
{
  double startTime = CPUTIME;
  HELP;

  argsdata args;
  args = parseArgs(argc, argv);

  nauty_check(WORDSIZE, MAXM, MAXN, NAUTYVERSIONID);
  initOptions();
  initCactiForests();

  calculate(&args);

  if (args.outputCount)
    printf(NSIZE_FORMAT " %lu %lu\n", args.n, a000080, a048781);
  if (!args.quiet)
    fprintf(stderr, ">Z a000080(%d)=%lu a048781(%d)=%lu in %3.2f sec\n", args.n, a000080, args.n, a048781, CPUTIME - startTime);

  return 0;
}
