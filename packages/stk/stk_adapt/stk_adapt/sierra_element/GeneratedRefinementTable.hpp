#ifndef STK_ADAPT_GENERATED_REFINEMENT_TABLES_HPP
#define STK_ADAPT_GENERATED_REFINEMENT_TABLES_HPP
/**  New ref topo info 
 *  ------------------
 *
 *  {Ord, Rnk-assoc, Ord-rnk-assoc, Ord-node-on-subcell, num-rnk-assoc, param-coord}
 *
 *   struct RefinementTopologyExtraEntry
 *   {
 *     unsigned ordinal_of_node;               // ordinal of node in the total list of nodes - corresponds to the shards node ordinal
 *     unsigned rank_of_subcell;               // rank of the subcell this node is associated with                                   
 *     unsigned ordinal_of_subcell;            // ordinal of the subcell in the shards numbering (e.g. edge # 3)
 *     unsigned ordinal_of_node_on_subcell;    // ordinal of the node on the subcell (whcih node it is on a subcell that has multiple nodes)
 *     unsigned num_nodes_on_subcell;          // how many nodes exist on the subcell                                                       
 *     double parametric_coordinates[3];
 *   };
 *       
 * Bootstrapping this file: to create this file, run the regression test RegressionTestUniformRefiner.cpp :: generate_tables after putting in
 *   a dummy entry in ./sierra_element/GeneratedRefinementTable.hpp.  The run will produce a local file, generated_refinement_tables.hpp 
 *   which can be checked against the gold copy of GeneratedRefinementTable.hpp, then copied over it.  Add a call below to generate the 
 *   actual new table data. 
 */



template<> RefTopoX RefinementTopologyExtra< shards:: Line<2>  > :: refinement_topology = {
  {	0,	0,	0,	0,	1,	{-1,	0,	0} },
  {	1,	0,	1,	0,	1,	{1,	0,	0} },
  {	2,	1,	0,	0,	1,	{0,	0,	0} } 

};

template<> RefTopoX RefinementTopologyExtra< shards:: ShellLine<2>  > :: refinement_topology = {
  {	0,	0,	0,	0,	1,	{-1,	0,	0} },
  {	1,	0,	1,	0,	1,	{1,	0,	0} },
  {	2,	1,	0,	0,	1,	{0,	0,	0} } 

};

template<> RefTopoX RefinementTopologyExtra< shards:: Quadrilateral<4>  > :: refinement_topology = {
  {	0,	0,	0,	0,	1,	{-1,	-1,	0} },
  {	1,	0,	1,	0,	1,	{1,	-1,	0} },
  {	2,	0,	2,	0,	1,	{1,	1,	0} },
  {	3,	0,	3,	0,	1,	{-1,	1,	0} },
  {	4,	1,	0,	0,	1,	{0,	-1,	0} },
  {	5,	1,	1,	0,	1,	{1,	0,	0} },
  {	6,	1,	2,	0,	1,	{0,	1,	0} },
  {	7,	1,	3,	0,	1,	{-1,	0,	0} },
  {	8,	2,	0,	0,	1,	{0,	0,	0} } 

};

template<> RefTopoX RefinementTopologyExtra< shards:: Triangle<3>  > :: refinement_topology = {
  {	0,	0,	0,	0,	1,	{0,	0,	0} },
  {	1,	0,	1,	0,	1,	{1,	0,	0} },
  {	2,	0,	2,	0,	1,	{0,	1,	0} },
  {	3,	1,	0,	0,	1,	{0.5,	0,	0} },
  {	4,	1,	1,	0,	1,	{0.5,	0.5,	0} },
  {	5,	1,	2,	0,	1,	{0,	0.5,	0} } 

};

template<> RefTopoX RefinementTopologyExtra< shards:: ShellTriangle<3>  > :: refinement_topology = {
  {	0,	0,	0,	0,	1,	{0,	0,	0} },
  {	1,	0,	1,	0,	1,	{1,	0,	0} },
  {	2,	0,	2,	0,	1,	{0,	1,	0} },
  {	3,	1,	0,	0,	1,	{0.5,	0,	0} },
  {	4,	1,	1,	0,	1,	{0.5,	0.5,	0} },
  {	5,	1,	2,	0,	1,	{0,	0.5,	0} } 

};

template<> RefTopoX RefinementTopologyExtra< shards:: ShellQuadrilateral<8>  > :: refinement_topology = {
  {	0,	0,	0,	0,	0,	{0,	0,	0} }

};

template<> RefTopoX RefinementTopologyExtra< shards:: ShellQuadrilateral<4>  > :: refinement_topology = {
  {	0,	0,	0,	0,	1,	{-1,	-1,	0} },
  {	1,	0,	1,	0,	1,	{1,	-1,	0} },
  {	2,	0,	2,	0,	1,	{1,	1,	0} },
  {	3,	0,	3,	0,	1,	{-1,	1,	0} },
  {	4,	1,	0,	0,	1,	{0,	-1,	0} },
  {	5,	1,	1,	0,	1,	{1,	0,	0} },
  {	6,	1,	2,	0,	1,	{0,	1,	0} },
  {	7,	1,	3,	0,	1,	{-1,	0,	0} },
  {	8,	2,	0,	0,	1,	{0,	0,	0} } 

};

template<> RefTopoX RefinementTopologyExtra< shards:: Tetrahedron<4>  > :: refinement_topology = {
  {	0,	0,	0,	0,	1,	{0,	0,	0} },
  {	1,	0,	1,	0,	1,	{1,	0,	0} },
  {	2,	0,	2,	0,	1,	{0,	1,	0} },
  {	3,	0,	3,	0,	1,	{0,	0,	1} },
  {	4,	1,	0,	0,	1,	{0.5,	0,	0} },
  {	5,	1,	1,	0,	1,	{0.5,	0.5,	0} },
  {	6,	1,	2,	0,	1,	{0,	0.5,	0} },
  {	7,	1,	3,	0,	1,	{0,	0,	0.5} },
  {	8,	1,	4,	0,	1,	{0.5,	0,	0.5} },
  {	9,	1,	5,	0,	1,	{0,	0.5,	0.5} } 

};

template<> RefTopoX RefinementTopologyExtra< shards:: Hexahedron<8>  > :: refinement_topology = {
  {	0,	0,	0,	0,	1,	{-1,	-1,	-1} },
  {	1,	0,	1,	0,	1,	{1,	-1,	-1} },
  {	2,	0,	2,	0,	1,	{1,	1,	-1} },
  {	3,	0,	3,	0,	1,	{-1,	1,	-1} },
  {	4,	0,	4,	0,	1,	{-1,	-1,	1} },
  {	5,	0,	5,	0,	1,	{1,	-1,	1} },
  {	6,	0,	6,	0,	1,	{1,	1,	1} },
  {	7,	0,	7,	0,	1,	{-1,	1,	1} },
  {	8,	1,	0,	0,	1,	{0,	-1,	-1} },
  {	9,	1,	1,	0,	1,	{1,	0,	-1} },
  {	10,	1,	2,	0,	1,	{0,	1,	-1} },
  {	11,	1,	3,	0,	1,	{-1,	0,	-1} },
  {	12,	1,	8,	0,	1,	{-1,	-1,	0} },
  {	13,	1,	9,	0,	1,	{1,	-1,	0} },
  {	14,	1,	10,	0,	1,	{1,	1,	0} },
  {	15,	1,	11,	0,	1,	{-1,	1,	0} },
  {	16,	1,	4,	0,	1,	{0,	-1,	1} },
  {	17,	1,	5,	0,	1,	{1,	0,	1} },
  {	18,	1,	6,	0,	1,	{0,	1,	1} },
  {	19,	1,	7,	0,	1,	{-1,	0,	1} },
  {	20,	3,	0,	0,	1,	{0,	0,	0} },
  {	21,	2,	4,	0,	1,	{0,	0,	-1} },
  {	22,	2,	5,	0,	1,	{0,	0,	1} },
  {	23,	2,	3,	0,	1,	{-1,	0,	0} },
  {	24,	2,	1,	0,	1,	{1,	0,	0} },
  {	25,	2,	0,	0,	1,	{0,	-1,	0} },
  {	26,	2,	2,	0,	1,	{0,	1,	0} } 

};

template<> RefTopoX RefinementTopologyExtra< shards:: Wedge<6>  > :: refinement_topology = {
  {	0,	0,	0,	0,	1,	{0,	0,	-1} },
  {	1,	0,	1,	0,	1,	{1,	0,	-1} },
  {	2,	0,	2,	0,	1,	{0,	1,	-1} },
  {	3,	0,	3,	0,	1,	{0,	0,	1} },
  {	4,	0,	4,	0,	1,	{1,	0,	1} },
  {	5,	0,	5,	0,	1,	{0,	1,	1} },
  {	6,	1,	0,	0,	1,	{0.5,	0,	-1} },
  {	7,	1,	1,	0,	1,	{0.5,	0.5,	-1} },
  {	8,	1,	2,	0,	1,	{0,	0.5,	-1} },
  {	9,	1,	6,	0,	1,	{0,	0,	0} },
  {	10,	1,	7,	0,	1,	{1,	0,	0} },
  {	11,	1,	8,	0,	1,	{0,	1,	0} },
  {	12,	1,	3,	0,	1,	{0.5,	0,	1} },
  {	13,	1,	4,	0,	1,	{0.5,	0.5,	1} },
  {	14,	1,	5,	0,	1,	{0,	0.5,	1} },
  {	15,	2,	0,	0,	1,	{0.5,	0,	0} },
  {	16,	2,	1,	0,	1,	{0.5,	0.5,	0} },
  {	17,	2,	2,	0,	1,	{0,	0.5,	0} } 

};

template<> RefTopoX RefinementTopologyExtra< shards:: Wedge<18>  > :: refinement_topology = {
  {	0,	0,	0,	0,	0,	{0,	0,	0} }

};

template<> RefTopoX RefinementTopologyExtra< shards:: Wedge<15>  > :: refinement_topology = {
  {	0,	0,	0,	0,	1,	{0,	0,	-1} },
  {	1,	0,	1,	0,	1,	{1,	0,	-1} },
  {	2,	0,	2,	0,	1,	{0,	1,	-1} },
  {	3,	0,	3,	0,	1,	{0,	0,	1} },
  {	4,	0,	4,	0,	1,	{1,	0,	1} },
  {	5,	0,	5,	0,	1,	{0,	1,	1} },
  {	6,	1,	0,	2,	3,	{0.5,	0,	-1} },
  {	7,	1,	1,	2,	3,	{0.5,	0.5,	-1} },
  {	8,	1,	2,	2,	3,	{0,	0.5,	-1} },
  {	9,	1,	6,	2,	3,	{0,	0,	0} },
  {	10,	1,	7,	2,	3,	{1,	0,	0} },
  {	11,	1,	8,	2,	3,	{0,	1,	0} },
  {	12,	1,	3,	2,	3,	{0.5,	0,	1} },
  {	13,	1,	4,	2,	3,	{0.5,	0.5,	1} },
  {	14,	1,	5,	2,	3,	{0,	0.5,	1} },
  {	15,	2,	0,	0,	5,	{0.5,	0,	0} },
  {	16,	2,	1,	0,	5,	{0.5,	0.5,	0} },
  {	17,	2,	2,	0,	5,	{0,	0.5,	0} },
  {	18,	1,	0,	0,	3,	{0.25,	0,	-1} },
  {	19,	1,	0,	1,	3,	{0.75,	0,	-1} },
  {	20,	1,	1,	0,	3,	{0.75,	0.25,	-1} },
  {	21,	1,	1,	1,	3,	{0.25,	0.75,	-1} },
  {	22,	1,	2,	0,	3,	{0,	0.75,	-1} },
  {	23,	1,	2,	1,	3,	{0,	0.25,	-1} },
  {	24,	2,	0,	4,	5,	{0.25,	0,	0} },
  {	25,	2,	0,	2,	5,	{0.75,	0,	0} },
  {	26,	2,	1,	4,	5,	{0.75,	0.25,	0} },
  {	27,	2,	1,	2,	5,	{0.25,	0.75,	0} },
  {	28,	2,	2,	3,	5,	{0,	0.75,	0} },
  {	29,	2,	2,	1,	5,	{0,	0.25,	0} },
  {	30,	1,	3,	0,	3,	{0.25,	0,	1} },
  {	31,	1,	3,	1,	3,	{0.75,	0,	1} },
  {	32,	1,	4,	0,	3,	{0.75,	0.25,	1} },
  {	33,	1,	4,	1,	3,	{0.25,	0.75,	1} },
  {	34,	1,	5,	0,	3,	{0,	0.75,	1} },
  {	35,	1,	5,	1,	3,	{0,	0.25,	1} },
  {	36,	2,	3,	1,	3,	{0.25,	0.5,	-1} },
  {	37,	2,	3,	0,	3,	{0.25,	0.25,	-1} },
  {	38,	2,	3,	2,	3,	{0.5,	0.25,	-1} },
  {	39,	3,	0,	0,	3,	{0.25,	0.5,	0} },
  {	40,	3,	0,	1,	3,	{0.25,	0.25,	0} },
  {	41,	3,	0,	2,	3,	{0.5,	0.25,	0} },
  {	42,	2,	4,	2,	3,	{0.25,	0.5,	1} },
  {	43,	2,	4,	0,	3,	{0.25,	0.25,	1} },
  {	44,	2,	4,	1,	3,	{0.5,	0.25,	1} },
  {	45,	1,	6,	0,	3,	{0,	0,	-0.5} },
  {	46,	2,	0,	1,	5,	{0.5,	0,	-0.5} },
  {	47,	1,	7,	0,	3,	{1,	0,	-0.5} },
  {	48,	2,	1,	1,	5,	{0.5,	0.5,	-0.5} },
  {	49,	1,	8,	0,	3,	{0,	1,	-0.5} },
  {	50,	2,	2,	4,	5,	{0,	0.5,	-0.5} },
  {	51,	1,	6,	1,	3,	{0,	0,	0.5} },
  {	52,	2,	0,	3,	5,	{0.5,	0,	0.5} },
  {	53,	1,	7,	1,	3,	{1,	0,	0.5} },
  {	54,	2,	1,	3,	5,	{0.5,	0.5,	0.5} },
  {	55,	1,	8,	1,	3,	{0,	1,	0.5} },
  {	56,	2,	2,	2,	5,	{0,	0.5,	0.5} } 

};

template<> RefTopoX RefinementTopologyExtra< shards:: Line<3>  > :: refinement_topology = {
  {	0,	0,	0,	0,	1,	{-1,	0,	0} },
  {	1,	0,	1,	0,	1,	{1,	0,	0} },
  {	2,	1,	0,	2,	3,	{0,	0,	0} },
  {	3,	1,	0,	0,	3,	{-0.5,	0,	0} },
  {	4,	1,	0,	1,	3,	{0.5,	0,	0} } 

};

template<> RefTopoX RefinementTopologyExtra< shards:: Triangle<6>  > :: refinement_topology = {
  {	0,	0,	0,	0,	1,	{0,	0,	0} },
  {	1,	0,	1,	0,	1,	{1,	0,	0} },
  {	2,	0,	2,	0,	1,	{0,	1,	0} },
  {	3,	1,	0,	2,	3,	{0.5,	0,	0} },
  {	4,	1,	1,	2,	3,	{0.5,	0.5,	0} },
  {	5,	1,	2,	2,	3,	{0,	0.5,	0} },
  {	6,	1,	0,	0,	3,	{0.25,	0,	0} },
  {	7,	1,	0,	1,	3,	{0.75,	0,	0} },
  {	8,	1,	1,	0,	3,	{0.75,	0.25,	0} },
  {	9,	1,	1,	1,	3,	{0.25,	0.75,	0} },
  {	10,	1,	2,	0,	3,	{0,	0.75,	0} },
  {	11,	1,	2,	1,	3,	{0,	0.25,	0} },
  {	12,	2,	0,	0,	3,	{0.25,	0.25,	0} },
  {	13,	2,	0,	1,	3,	{0.5,	0.25,	0} },
  {	14,	2,	0,	2,	3,	{0.25,	0.5,	0} } 

};

template<> RefTopoX RefinementTopologyExtra< shards:: Quadrilateral<8>  > :: refinement_topology = {
  {	0,	0,	0,	0,	1,	{-1,	-1,	0} },
  {	1,	0,	1,	0,	1,	{1,	-1,	0} },
  {	2,	0,	2,	0,	1,	{1,	1,	0} },
  {	3,	0,	3,	0,	1,	{-1,	1,	0} },
  {	4,	1,	0,	2,	3,	{0,	-1,	0} },
  {	5,	1,	1,	2,	3,	{1,	0,	0} },
  {	6,	1,	2,	2,	3,	{0,	1,	0} },
  {	7,	1,	3,	2,	3,	{-1,	0,	0} },
  {	8,	2,	0,	8,	9,	{0,	0,	0} },
  {	9,	1,	0,	0,	3,	{-0.5,	-1,	0} },
  {	10,	1,	0,	1,	3,	{0.5,	-1,	0} },
  {	11,	1,	1,	0,	3,	{1,	-0.5,	0} },
  {	12,	1,	1,	1,	3,	{1,	0.5,	0} },
  {	13,	1,	2,	0,	3,	{0.5,	1,	0} },
  {	14,	1,	2,	1,	3,	{-0.5,	1,	0} },
  {	15,	1,	3,	0,	3,	{-1,	0.5,	0} },
  {	16,	1,	3,	1,	3,	{-1,	-0.5,	0} },
  {	17,	2,	0,	4,	9,	{0,	-0.5,	0} },
  {	18,	2,	0,	5,	9,	{0.5,	0,	0} },
  {	19,	2,	0,	6,	9,	{0,	0.5,	0} },
  {	20,	2,	0,	7,	9,	{-0.5,	0,	0} } 

};

template<> RefTopoX RefinementTopologyExtra< shards:: Quadrilateral<9>  > :: refinement_topology = {
  {	0,	0,	0,	0,	1,	{-1,	-1,	0} },
  {	1,	0,	1,	0,	1,	{1,	-1,	0} },
  {	2,	0,	2,	0,	1,	{1,	1,	0} },
  {	3,	0,	3,	0,	1,	{-1,	1,	0} },
  {	4,	1,	0,	2,	3,	{0,	-1,	0} },
  {	5,	1,	1,	2,	3,	{1,	0,	0} },
  {	6,	1,	2,	2,	3,	{0,	1,	0} },
  {	7,	1,	3,	2,	3,	{-1,	0,	0} },
  {	8,	2,	0,	8,	9,	{0,	0,	0} },
  {	9,	1,	0,	0,	3,	{-0.5,	-1,	0} },
  {	10,	1,	0,	1,	3,	{0.5,	-1,	0} },
  {	11,	1,	1,	0,	3,	{1,	-0.5,	0} },
  {	12,	1,	1,	1,	3,	{1,	0.5,	0} },
  {	13,	1,	2,	0,	3,	{0.5,	1,	0} },
  {	14,	1,	2,	1,	3,	{-0.5,	1,	0} },
  {	15,	1,	3,	0,	3,	{-1,	0.5,	0} },
  {	16,	1,	3,	1,	3,	{-1,	-0.5,	0} },
  {	17,	2,	0,	4,	9,	{0,	-0.5,	0} },
  {	18,	2,	0,	5,	9,	{0.5,	0,	0} },
  {	19,	2,	0,	6,	9,	{0,	0.5,	0} },
  {	20,	2,	0,	7,	9,	{-0.5,	0,	0} },
  {	21,	2,	0,	0,	9,	{-0.5,	-0.5,	0} },
  {	22,	2,	0,	1,	9,	{0.5,	-0.5,	0} },
  {	23,	2,	0,	2,	9,	{0.5,	0.5,	0} },
  {	24,	2,	0,	3,	9,	{-0.5,	0.5,	0} } 

};


template<> RefTopoX RefinementTopologyExtra< shards:: Hexahedron<27>  > :: refinement_topology = {
  {	0,	0,	0,	0,	1,	{-1,	-1,	-1} },
  {	1,	0,	1,	0,	1,	{1,	-1,	-1} },
  {	2,	0,	2,	0,	1,	{1,	1,	-1} },
  {	3,	0,	3,	0,	1,	{-1,	1,	-1} },
  {	4,	0,	4,	0,	1,	{-1,	-1,	1} },
  {	5,	0,	5,	0,	1,	{1,	-1,	1} },
  {	6,	0,	6,	0,	1,	{1,	1,	1} },
  {	7,	0,	7,	0,	1,	{-1,	1,	1} },
  {	8,	1,	0,	2,	3,	{0,	-1,	-1} },
  {	9,	1,	1,	2,	3,	{1,	0,	-1} },
  {	10,	1,	2,	2,	3,	{0,	1,	-1} },
  {	11,	1,	3,	2,	3,	{-1,	0,	-1} },
  {	12,	1,	8,	2,	3,	{-1,	-1,	0} },
  {	13,	1,	9,	2,	3,	{1,	-1,	0} },
  {	14,	1,	10,	2,	3,	{1,	1,	0} },
  {	15,	1,	11,	2,	3,	{-1,	1,	0} },
  {	16,	1,	4,	2,	3,	{0,	-1,	1} },
  {	17,	1,	5,	2,	3,	{1,	0,	1} },
  {	18,	1,	6,	2,	3,	{0,	1,	1} },
  {	19,	1,	7,	2,	3,	{-1,	0,	1} },
  {	20,	3,	0,	0,	27,	{0,	0,	0} },
  {	21,	2,	4,	8,	9,	{0,	0,	-1} },
  {	22,	2,	5,	8,	9,	{0,	0,	1} },
  {	23,	2,	3,	8,	9,	{-1,	0,	0} },
  {	24,	2,	1,	8,	9,	{1,	0,	0} },
  {	25,	2,	0,	8,	9,	{0,	-1,	0} },
  {	26,	2,	2,	8,	9,	{0,	1,	0} },
  {	27,	1,	0,	0,	3,	{-0.5,	-1,	-1} },
  {	28,	1,	0,	1,	3,	{0.5,	-1,	-1} },
  {	29,	1,	1,	0,	3,	{1,	-0.5,	-1} },
  {	30,	1,	1,	1,	3,	{1,	0.5,	-1} },
  {	31,	1,	2,	0,	3,	{0.5,	1,	-1} },
  {	32,	1,	2,	1,	3,	{-0.5,	1,	-1} },
  {	33,	1,	3,	0,	3,	{-1,	0.5,	-1} },
  {	34,	1,	3,	1,	3,	{-1,	-0.5,	-1} },
  {	35,	1,	8,	0,	3,	{-1,	-1,	-0.5} },
  {	36,	1,	9,	0,	3,	{1,	-1,	-0.5} },
  {	37,	1,	10,	0,	3,	{1,	1,	-0.5} },
  {	38,	1,	11,	0,	3,	{-1,	1,	-0.5} },
  {	39,	1,	8,	1,	3,	{-1,	-1,	0.5} },
  {	40,	1,	9,	1,	3,	{1,	-1,	0.5} },
  {	41,	1,	10,	1,	3,	{1,	1,	0.5} },
  {	42,	1,	11,	1,	3,	{-1,	1,	0.5} },
  {	43,	1,	4,	0,	3,	{-0.5,	-1,	1} },
  {	44,	1,	4,	1,	3,	{0.5,	-1,	1} },
  {	45,	1,	5,	0,	3,	{1,	-0.5,	1} },
  {	46,	1,	5,	1,	3,	{1,	0.5,	1} },
  {	47,	1,	6,	0,	3,	{0.5,	1,	1} },
  {	48,	1,	6,	1,	3,	{-0.5,	1,	1} },
  {	49,	1,	7,	0,	3,	{-1,	0.5,	1} },
  {	50,	1,	7,	1,	3,	{-1,	-0.5,	1} },
  {	51,	2,	0,	7,	9,	{-0.5,	-1,	0} },
  {	52,	2,	0,	5,	9,	{0.5,	-1,	0} },
  {	53,	2,	1,	7,	9,	{1,	-0.5,	0} },
  {	54,	2,	1,	5,	9,	{1,	0.5,	0} },
  {	55,	2,	2,	7,	9,	{0.5,	1,	0} },
  {	56,	2,	2,	5,	9,	{-0.5,	1,	0} },
  {	57,	2,	3,	6,	9,	{-1,	0.5,	0} },
  {	58,	2,	3,	4,	9,	{-1,	-0.5,	0} },
  {	59,	2,	0,	4,	9,	{0,	-1,	-0.5} },
  {	60,	2,	4,	7,	9,	{0,	-0.5,	-1} },
  {	61,	2,	4,	5,	9,	{0,	0.5,	-1} },
  {	62,	2,	2,	4,	9,	{0,	1,	-0.5} },
  {	63,	2,	2,	6,	9,	{0,	1,	0.5} },
  {	64,	2,	5,	6,	9,	{0,	0.5,	1} },
  {	65,	2,	5,	4,	9,	{0,	-0.5,	1} },
  {	66,	2,	0,	6,	9,	{0,	-1,	0.5} },
  {	67,	2,	4,	4,	9,	{-0.5,	0,	-1} },
  {	68,	2,	4,	6,	9,	{0.5,	0,	-1} },
  {	69,	2,	1,	4,	9,	{1,	0,	-0.5} },
  {	70,	2,	1,	6,	9,	{1,	0,	0.5} },
  {	71,	2,	5,	5,	9,	{0.5,	0,	1} },
  {	72,	2,	5,	7,	9,	{-0.5,	0,	1} },
  {	73,	2,	3,	5,	9,	{-1,	0,	0.5} },
  {	74,	2,	3,	7,	9,	{-1,	0,	-0.5} },
  {	75,	3,	0,	1,	27,	{0,	-0.5,	0} },
  {	76,	3,	0,	2,	27,	{0,	0.5,	0} },
  {	77,	3,	0,	3,	27,	{-0.5,	0,	0} },
  {	78,	3,	0,	4,	27,	{0.5,	0,	0} },
  {	79,	3,	0,	5,	27,	{0,	0,	-0.5} },
  {	80,	3,	0,	6,	27,	{0,	0,	0.5} },
  {	81,	3,	0,	7,	27,	{0,	0,	0} },
  {	82,	3,	0,	8,	27,	{0,	0,	0} },
  {	83,	3,	0,	9,	27,	{0,	0,	0} },
  {	84,	3,	0,	10,	27,	{0,	0,	0} },
  {	85,	3,	0,	11,	27,	{0,	0,	0} },
  {	86,	3,	0,	12,	27,	{0,	0,	0} },
  {	87,	3,	0,	13,	27,	{0,	0,	0} },
  {	88,	3,	0,	14,	27,	{0,	0,	0} },
  {	89,	2,	4,	0,	9,	{-0.5,	-0.5,	-1} },
  {	90,	2,	4,	1,	9,	{-0.5,	0.5,	-1} },
  {	91,	2,	4,	2,	9,	{0.5,	0.5,	-1} },
  {	92,	2,	4,	3,	9,	{0.5,	-0.5,	-1} },
  {	93,	2,	5,	0,	9,	{-0.5,	-0.5,	1} },
  {	94,	2,	5,	1,	9,	{0.5,	-0.5,	1} },
  {	95,	2,	5,	2,	9,	{0.5,	0.5,	1} },
  {	96,	2,	5,	3,	9,	{-0.5,	0.5,	1} },
  {	97,	2,	3,	0,	9,	{-1,	-0.5,	-0.5} },
  {	98,	2,	3,	1,	9,	{-1,	-0.5,	0.5} },
  {	99,	2,	3,	2,	9,	{-1,	0.5,	0.5} },
  {	100,	2,	3,	3,	9,	{-1,	0.5,	-0.5} },
  {	101,	2,	1,	0,	9,	{1,	-0.5,	-0.5} },
  {	102,	2,	1,	1,	9,	{1,	0.5,	-0.5} },
  {	103,	2,	1,	2,	9,	{1,	0.5,	0.5} },
  {	104,	2,	1,	3,	9,	{1,	-0.5,	0.5} },
  {	105,	2,	0,	0,	9,	{-0.5,	-1,	-0.5} },
  {	106,	2,	0,	1,	9,	{0.5,	-1,	-0.5} },
  {	107,	2,	0,	2,	9,	{0.5,	-1,	0.5} },
  {	108,	2,	0,	3,	9,	{-0.5,	-1,	0.5} },
  {	109,	2,	2,	0,	9,	{0.5,	1,	-0.5} },
  {	110,	2,	2,	1,	9,	{-0.5,	1,	-0.5} },
  {	111,	2,	2,	2,	9,	{-0.5,	1,	0.5} },
  {	112,	2,	2,	3,	9,	{0.5,	1,	0.5} },
  {	113,	3,	0,	15,	27,	{0,	-0.5,	-0.5} },
  {	114,	3,	0,	16,	27,	{0,	0.5,	-0.5} },
  {	115,	3,	0,	17,	27,	{0,	0.5,	0.5} },
  {	116,	3,	0,	18,	27,	{0,	-0.5,	0.5} },
  {	117,	3,	0,	19,	27,	{-0.5,	-0.5,	0} },
  {	118,	3,	0,	20,	27,	{0.5,	-0.5,	0} },
  {	119,	3,	0,	21,	27,	{0.5,	0.5,	0} },
  {	120,	3,	0,	22,	27,	{-0.5,	0.5,	0} },
  {	121,	3,	0,	23,	27,	{-0.5,	0,	-0.5} },
  {	122,	3,	0,	24,	27,	{0.5,	0,	-0.5} },
  {	123,	3,	0,	25,	27,	{0.5,	0,	0.5} },
  {	124,	3,	0,	26,	27,	{-0.5,	0,	0.5} } 

};

template<> RefTopoX RefinementTopologyExtra< shards:: Hexahedron<20>  > :: refinement_topology = {
  {	0,	0,	0,	0,	1,	{-1,	-1,	-1} },
  {	1,	0,	1,	0,	1,	{1,	-1,	-1} },
  {	2,	0,	2,	0,	1,	{1,	1,	-1} },
  {	3,	0,	3,	0,	1,	{-1,	1,	-1} },
  {	4,	0,	4,	0,	1,	{-1,	-1,	1} },
  {	5,	0,	5,	0,	1,	{1,	-1,	1} },
  {	6,	0,	6,	0,	1,	{1,	1,	1} },
  {	7,	0,	7,	0,	1,	{-1,	1,	1} },
  {	8,	1,	0,	2,	3,	{0,	-1,	-1} },
  {	9,	1,	1,	2,	3,	{1,	0,	-1} },
  {	10,	1,	2,	2,	3,	{0,	1,	-1} },
  {	11,	1,	3,	2,	3,	{-1,	0,	-1} },
  {	12,	1,	8,	2,	3,	{-1,	-1,	0} },
  {	13,	1,	9,	2,	3,	{1,	-1,	0} },
  {	14,	1,	10,	2,	3,	{1,	1,	0} },
  {	15,	1,	11,	2,	3,	{-1,	1,	0} },
  {	16,	1,	4,	2,	3,	{0,	-1,	1} },
  {	17,	1,	5,	2,	3,	{1,	0,	1} },
  {	18,	1,	6,	2,	3,	{0,	1,	1} },
  {	19,	1,	7,	2,	3,	{-1,	0,	1} },
  {	20,	3,	0,	0,	7,	{0,	0,	0} },
  {	21,	2,	4,	8,	9,	{0,	0,	-1} },
  {	22,	2,	5,	8,	9,	{0,	0,	1} },
  {	23,	2,	3,	8,	9,	{-1,	0,	0} },
  {	24,	2,	1,	8,	9,	{1,	0,	0} },
  {	25,	2,	0,	8,	9,	{0,	-1,	0} },
  {	26,	2,	2,	8,	9,	{0,	1,	0} },
  {	27,	1,	0,	0,	3,	{-0.5,	-1,	-1} },
  {	28,	1,	0,	1,	3,	{0.5,	-1,	-1} },
  {	29,	1,	1,	0,	3,	{1,	-0.5,	-1} },
  {	30,	1,	1,	1,	3,	{1,	0.5,	-1} },
  {	31,	1,	2,	0,	3,	{0.5,	1,	-1} },
  {	32,	1,	2,	1,	3,	{-0.5,	1,	-1} },
  {	33,	1,	3,	0,	3,	{-1,	0.5,	-1} },
  {	34,	1,	3,	1,	3,	{-1,	-0.5,	-1} },
  {	35,	1,	8,	0,	3,	{-1,	-1,	-0.5} },
  {	36,	1,	9,	0,	3,	{1,	-1,	-0.5} },
  {	37,	1,	10,	0,	3,	{1,	1,	-0.5} },
  {	38,	1,	11,	0,	3,	{-1,	1,	-0.5} },
  {	39,	1,	8,	1,	3,	{-1,	-1,	0.5} },
  {	40,	1,	9,	1,	3,	{1,	-1,	0.5} },
  {	41,	1,	10,	1,	3,	{1,	1,	0.5} },
  {	42,	1,	11,	1,	3,	{-1,	1,	0.5} },
  {	43,	1,	4,	0,	3,	{-0.5,	-1,	1} },
  {	44,	1,	4,	1,	3,	{0.5,	-1,	1} },
  {	45,	1,	5,	0,	3,	{1,	-0.5,	1} },
  {	46,	1,	5,	1,	3,	{1,	0.5,	1} },
  {	47,	1,	6,	0,	3,	{0.5,	1,	1} },
  {	48,	1,	6,	1,	3,	{-0.5,	1,	1} },
  {	49,	1,	7,	0,	3,	{-1,	0.5,	1} },
  {	50,	1,	7,	1,	3,	{-1,	-0.5,	1} },
  {	51,	2,	0,	7,	9,	{-0.5,	-1,	0} },
  {	52,	2,	0,	5,	9,	{0.5,	-1,	0} },
  {	53,	2,	1,	7,	9,	{1,	-0.5,	0} },
  {	54,	2,	1,	5,	9,	{1,	0.5,	0} },
  {	55,	2,	2,	7,	9,	{0.5,	1,	0} },
  {	56,	2,	2,	5,	9,	{-0.5,	1,	0} },
  {	57,	2,	3,	6,	9,	{-1,	0.5,	0} },
  {	58,	2,	3,	4,	9,	{-1,	-0.5,	0} },
  {	59,	2,	0,	4,	9,	{0,	-1,	-0.5} },
  {	60,	2,	4,	7,	9,	{0,	-0.5,	-1} },
  {	61,	2,	4,	5,	9,	{0,	0.5,	-1} },
  {	62,	2,	2,	4,	9,	{0,	1,	-0.5} },
  {	63,	2,	2,	6,	9,	{0,	1,	0.5} },
  {	64,	2,	5,	6,	9,	{0,	0.5,	1} },
  {	65,	2,	5,	4,	9,	{0,	-0.5,	1} },
  {	66,	2,	0,	6,	9,	{0,	-1,	0.5} },
  {	67,	2,	4,	4,	9,	{-0.5,	0,	-1} },
  {	68,	2,	4,	6,	9,	{0.5,	0,	-1} },
  {	69,	2,	1,	4,	9,	{1,	0,	-0.5} },
  {	70,	2,	1,	6,	9,	{1,	0,	0.5} },
  {	71,	2,	5,	5,	9,	{0.5,	0,	1} },
  {	72,	2,	5,	7,	9,	{-0.5,	0,	1} },
  {	73,	2,	3,	5,	9,	{-1,	0,	0.5} },
  {	74,	2,	3,	7,	9,	{-1,	0,	-0.5} },
  {	75,	3,	0,	1,	7,	{0,	-0.5,	0} },
  {	76,	3,	0,	2,	7,	{0,	0.5,	0} },
  {	77,	3,	0,	3,	7,	{-0.5,	0,	0} },
  {	78,	3,	0,	4,	7,	{0.5,	0,	0} },
  {	79,	3,	0,	5,	7,	{0,	0,	-0.5} },
  {	80,	3,	0,	6,	7,	{0,	0,	0.5} } 

};

template<> RefTopoX RefinementTopologyExtra< shards:: Tetrahedron<10>  > :: refinement_topology = {
  {	0,	0,	0,	0,	1,	{0,	0,	0} },
  {	1,	0,	1,	0,	1,	{1,	0,	0} },
  {	2,	0,	2,	0,	1,	{0,	1,	0} },
  {	3,	0,	3,	0,	1,	{0,	0,	1} },
  {	4,	1,	0,	2,	3,	{0.5,	0,	0} },
  {	5,	1,	1,	2,	3,	{0.5,	0.5,	0} },
  {	6,	1,	2,	2,	3,	{0,	0.5,	0} },
  {	7,	1,	3,	2,	3,	{0,	0,	0.5} },
  {	8,	1,	4,	2,	3,	{0.5,	0,	0.5} },
  {	9,	1,	5,	2,	3,	{0,	0.5,	0.5} },
  {	10,	1,	0,	0,	3,	{0.25,	0,	0} },
  {	11,	1,	0,	1,	3,	{0.75,	0,	0} },
  {	12,	1,	1,	0,	3,	{0.75,	0.25,	0} },
  {	13,	1,	1,	1,	3,	{0.25,	0.75,	0} },
  {	14,	1,	2,	0,	3,	{0,	0.75,	0} },
  {	15,	1,	2,	1,	3,	{0,	0.25,	0} },
  {	16,	1,	3,	0,	3,	{0,	0,	0.25} },
  {	17,	1,	4,	0,	3,	{0.75,	0,	0.25} },
  {	18,	1,	5,	0,	3,	{0,	0.75,	0.25} },
  {	19,	1,	3,	1,	3,	{0,	0,	0.75} },
  {	20,	1,	4,	1,	3,	{0.25,	0,	0.75} },
  {	21,	1,	5,	1,	3,	{0,	0.25,	0.75} },
  {	22,	3,	0,	0,	1,	{0.25,	0.25,	0.25} },
  {	23,	2,	0,	2,	3,	{0.25,	0,	0.5} },
  {	24,	2,	3,	1,	3,	{0.25,	0.5,	0} },
  {	25,	2,	2,	0,	3,	{0,	0.25,	0.25} },
  {	26,	2,	1,	0,	3,	{0.5,	0.25,	0.25} },
  {	27,	2,	0,	0,	3,	{0.25,	0,	0.25} },
  {	28,	2,	1,	1,	3,	{0.25,	0.5,	0.25} },
  {	29,	2,	2,	1,	3,	{0,	0.25,	0.5} },
  {	30,	2,	3,	2,	3,	{0.5,	0.25,	0} },
  {	31,	2,	3,	0,	3,	{0.25,	0.25,	0} },
  {	32,	2,	1,	2,	3,	{0.25,	0.25,	0.5} },
  {	33,	2,	2,	2,	3,	{0,	0.5,	0.25} },
  {	34,	2,	0,	1,	3,	{0.5,	0,	0.25} } 

};
#endif
