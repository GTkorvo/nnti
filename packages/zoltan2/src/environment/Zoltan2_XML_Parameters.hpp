// @HEADER
//
// ***********************************************************************
//
//   Zoltan2: A package of combinatorial algorithms for scientific computing
//                  Copyright 2012 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Karen Devine      (kddevin@sandia.gov)
//                    Erik Boman        (egboman@sandia.gov)
//                    Siva Rajamanickam (srajama@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
// 
// This file was automatically generated by CMake
// with the following command:
// ./xmlToHeaderDefinition /Users/kddevin/code/Trilinos/packages/zoltan2/data/parameters.xml /Users/kddevin/code/Trilinos/Obj_zoltan2/packages/zoltan2/src/Zoltan2_XML_Parameters.hpp.
// 
#ifndef ZOLTAN2_PARAMETER_DEFINITION_HEADER
#define ZOLTAN2_PARAMETER_DEFINITION_HEADER

#define ZOLTAN2_XML_PARAMETER_STRING " \
  <ParameterList name=\"zoltan2ValidatingParameters\"> \
   <Parameter  \
    id=\"0\" name=\"error_check_level\" type=\"string\" validatorId=\"0\" value=\"basic_assertions\" \
    docString='  the amount of error checking performed \
    (If the compile flag Z2_OMIT_ALL_ERROR_CHECKING was set, \
    then error checking code is not executed at runtime.)' \
    /> \
   <Parameter  \
    id=\"1\" name=\"debug_level\" type=\"string\" validatorId=\"1\" value=\"basic_status\" \
    docString='  the amount of status/debugging/warning information to print' \
    /> \
   <Parameter  \
    id=\"2\" name=\"timer_type\" type=\"string\" validatorId=\"2\" value=\"no_timers\" \
    docString='  the type of timing information to collect \
    (If the compile flag Z2_OMIT_ALL_PROFILING was set, \
    then the timing code is not executed at runtime.)' \
    /> \
   <Parameter  \
    id=\"3\" name=\"debug_output_stream\" type=\"string\" validatorId=\"3\" value=\"cout\" \
    docString='  output stream for debug/status/warning messages (default cout)' \
    /> \
   <Parameter  \
    id=\"4\" name=\"timer_output_stream\" type=\"string\" validatorId=\"4\" value=\"cout\" \
    docString='  output stream for timing report (default cout)' \
    /> \
   <Parameter  \
    id=\"5\" name=\"memory_output_stream\" type=\"string\" validatorId=\"5\" value=\"cout\" \
    docString='  output stream for memory usage messages (default cout)' \
    /> \
   <Parameter  \
    id=\"6\" name=\"debug_output_file\" type=\"string\" validatorId=\"6\" value=\"/dev/null\" \
    docString='  name of file to which debug/status messages should be written \
    (process rank will be included in file name)' \
    /> \
   <Parameter  \
    id=\"7\" name=\"timer_output_file\" type=\"string\" validatorId=\"6\" value=\"/dev/null\" \
    docString='  name of file to which timing information should be written \
    (process rank will be included in file name)' \
    /> \
   <Parameter  \
    id=\"8\" name=\"memory_output_file\" type=\"string\" validatorId=\"6\" value=\"/dev/null\" \
    docString='  name of file to which memory profiling information should be written \
    (process rank will be included in file name)' \
    /> \
   <Parameter  \
    id=\"9\" name=\"debug_procs\" type=\"string\" validatorId=\"7\" value=\"0\" \
    docString='  list of ranks that output debugging/status messages (default \"0\")' \
    /> \
   <Parameter  \
    id=\"10\" name=\"pqParts\" type=\"string\" validatorId=\"8\" value=\"0\" \
    docString='  list of parts for multiJagged partitioning algorithm. As many as the dimension count.' \
    /> \
   <Parameter  \
    id=\"11\" name=\"memory_procs\" type=\"string\" validatorId=\"9\" value=\"0\" \
    docString='  list of ranks that do memory profiling information (default \"0\")' \
    /> \
   <Parameter  \
    id=\"12\" name=\"speed_versus_quality\" type=\"string\" validatorId=\"10\" value=\"balance\" \
    docString='  When algorithm choices exist, opt for speed or solution quality? \
    (Default is a balance of speed and quality)' \
    /> \
   <Parameter  \
    id=\"13\" name=\"memory_versus_speed\" type=\"string\" validatorId=\"11\" value=\"balance\" \
    docString='  When algorithm choices exist, opt for the use of less memory \
    at the expense of runtime \
    (Default is a balance of memory conservation and speed)' \
    /> \
   <Parameter  \
    id=\"14\" name=\"random_seed\" type=\"string\" validatorId=\"12\" value=\"0.5\" \
    docString='  random seed' \
    /> \
   <Parameter  \
    id=\"15\" name=\"order_method\" type=\"string\" validatorId=\"13\" value=\"rcm\" \
    docString='  order algorithm' \
    /> \
   <Parameter  \
    id=\"16\" name=\"order_package\" type=\"string\" validatorId=\"14\" value=\"amd\" \
    docString='  package to use in ordering' \
    /> \
   <Parameter  \
    id=\"17\" name=\"compute_metrics\" type=\"string\" validatorId=\"15\" value=\"no\" \
    docString='  Compute metrics after computing solution' \
    /> \
   <Parameter  \
    id=\"18\" name=\"topology\" type=\"string\" validatorId=\"16\" value=\"\" \
    docString='Topology of node to be used in hierarchical partitioning \
    \"2,4\"  for dual-socket quad-core \
    \"2,2,6\"  for dual-socket, dual-die, six-core \
    \"2,2,3\"  for dual-socket, dual-die, six-core but \
               with only three partitions per die' \
    /> \
   <Parameter  \
    id=\"19\" name=\"randomize_input\" type=\"string\" validatorId=\"17\" value=\"no\" \
    docString='  randomize input prior to partitioning' \
    /> \
   <Parameter  \
    id=\"20\" name=\"partitioning_objective\" type=\"string\" validatorId=\"18\" value=\"balance_object_weight\" \
    docString='  objective of partitioning (default depends on algorithm)' \
    /> \
   <Parameter  \
    id=\"21\" name=\"imbalance_tolerance\" type=\"string\" validatorId=\"19\" value=\"1.1\" \
    docString='  imbalance tolerance, ratio of maximum load over average load (default 1.1)' \
    /> \
   <Parameter  \
    id=\"22\" name=\"num_global_parts\" type=\"string\" validatorId=\"20\" value=\"0\" \
    docString='  global number of parts to compute (default is number of processes)' \
    /> \
   <Parameter  \
    id=\"23\" name=\"num_local_parts\" type=\"string\" validatorId=\"21\" value=\"0\" \
    docString='  number of parts to compute for this process(default is one)' \
    /> \
   <Parameter  \
    id=\"24\" name=\"partitioning_approach\" type=\"string\" validatorId=\"22\" value=\"partition\" \
    docString='  Partition from scratch, partition incrementally from current \
    partition, of partition from scratch but maximize overlap \
    with the current partition (default is \"partition\" from scratch)' \
    /> \
   <Parameter  \
    id=\"25\" name=\"objects_to_partition\" type=\"string\" validatorId=\"23\" value=\"graph_vertices\" \
    docString='  Objects to be partitioned (defaults are \"matrix_rows\" for \
    matrix input, \"mesh_nodes\" for mesh input, and \"graph_vertices\" \
    for graph input)' \
    /> \
   <Parameter  \
    id=\"26\" name=\"model\" type=\"string\" validatorId=\"24\" value=\"graph\" \
    docString='  This is a low level parameter.  Normally the library will choose \
    a computational model based on the algorithm or objective specified \
    by the user.' \
    /> \
   <Parameter  \
    id=\"27\" name=\"algorithm\" type=\"string\" validatorId=\"25\" value=\"random\" \
    docString='  partitioning algorithm' \
    /> \
   <Parameter  \
    id=\"28\" name=\"rectilinear_blocks\" type=\"string\" validatorId=\"26\" value=\"no\" \
    docString='  If true, then when a cut is made, all of the dots located on the cut \
    are moved to the same side of the cut. The resulting regions are then \
    rectilinear.  The resulting load balance may not be as good as when \
    the group of dots is split by the cut. Default is false.' \
    /> \
   <Parameter  \
    id=\"29\" name=\"average_cuts\" type=\"string\" validatorId=\"27\" value=\"no\" \
    docString='  When true, coordinates of RCB cutting planes are computed to be  \
    the average of the coordinates of the closest object on each side  \
    of the cut. Otherwise, coordinates of cutting planes may equal  \
    those of one of the closest objects. Default is false.' \
    /> \
   <Parameter  \
    id=\"30\" name=\"bisection_num_test_cuts\" type=\"int\" validatorId=\"28\" value=\"1\" \
    docString='  Experimental: number of test cuts to do simultaneously (default is 1)' \
    /> \
   <Parameter  \
    id=\"31\" name=\"symmetrize_input\" type=\"string\" validatorId=\"29\" value=\"no\" \
    docString='  Symmetrize input prior to pList.  If \"transpose\", \
    symmetrize A by computing A plus ATranspose.  If \"bipartite\", \
    A becomes [[0 A][ATranspose 0]].  ' \
    /> \
   <Parameter  \
    id=\"32\" name=\"subset_graph\" type=\"string\" validatorId=\"30\" value=\"no\" \
    docString='  If \"yes\", the graph input is to be subsetted.  If a vertex neighbor \
    is not a valid vertex, it will be omitted from the pList.  Otherwise, \
    an invalid neighbor identifier is considered an error.' \
    /> \
   <Parameter  \
    id=\"35\" name=\"parallel_part_calculation_count\" type=\"int\" validatorId=\"33\" value=\"1\" \
    docString=\"The number of parts whose cut coordinates will be calculated concurently.\" \
    /> \
   <Parameter  \
    id=\"36\" name=\"migration_imbalance_cut_off\" type=\"string\" validatorId=\"19\" value=\"1.1\" \
    docString='  migration_imbalance_cut_off, the minimum imbalance of the processors to avoid migration (default 1.1)' \
    /> \
   <Parameter  \
    id=\"37\" name=\"migration_all_to_all_type\" type=\"int\" validatorId=\"37\" value=\"1\" \
    docString=\"Migration type, 0 for naive migration, 1 for smarter migration.\" \
    /> \
   <Parameter  \
    id=\"38\" name=\"migration_check_option\" type=\"int\" validatorId=\"38\" value=\"1\" \
    docString=\"Migration option, 0 for decision depending on the imbalance, 1 for forcing migration, 2 for avoiding migration\" \
    /> \
   <Parameter  \
    id=\"39\" name=\"migration_processor_assignment_type\" type=\"int\" validatorId=\"39\" value=\"1\" \
    docString=\"Migration processor assignment type, 0 for assignning procs with respect to weight, otherwise for assigning procs with respect to closeness\" \
    /> \
   <Parameter  \
    id=\"40\" name=\"remap_parts\" type=\"string\" validatorId=\"40\" value=\"no\" \
    docString='  remap part numbers to minimize migration between old and new partitions' \
    /> \
    <Validators> \
      <Validator defaultParameterName=\"error_check_level\" integralValue=\"int\" type=\"StringIntegralValidator(int)\" validatorId=\"0\"> \
        <String integralValue=\"0\" stringDoc=\"no assertions will be performed\" stringValue=\"no_assertions\"/> \
        <String integralValue=\"1\" stringDoc=\"typical checks of argument validity (fast, default)\" stringValue=\"basic_assertions\"/> \
        <String integralValue=\"2\" stringDoc=\"additional checks, i.e. is input graph a valid graph)\" stringValue=\"complex_assertions\"/> \
        <String integralValue=\"3\" stringDoc=\"check for everything including logic errors (slowest)\" stringValue=\"debug_mode_assertions\"/> \
      </Validator> \
      <Validator defaultParameterName=\"debug_output_stream\" integralValue=\"int\" type=\"StringIntegralValidator(int)\" validatorId=\"3\"> \
        <String integralValue=\"0\" stringValue=\"std::cout\"/> \
        <String integralValue=\"0\" stringValue=\"cout\"/> \
        <String integralValue=\"0\" stringValue=\"stdout\"/> \
        <String integralValue=\"1\" stringValue=\"std::cerr\"/> \
        <String integralValue=\"1\" stringValue=\"cerr\"/> \
        <String integralValue=\"1\" stringValue=\"stderr\"/> \
        <String integralValue=\"2\" stringValue=\"/dev/null\"/> \
        <String integralValue=\"2\" stringValue=\"null\"/> \
      </Validator> \
      <Validator defaultParameterName=\"debug_level\" integralValue=\"int\" type=\"StringIntegralValidator(int)\" validatorId=\"1\"> \
        <String integralValue=\"0\" stringDoc=\"library outputs no status information\" stringValue=\"no_status\"/> \
        <String integralValue=\"1\" stringDoc=\"library outputs basic status information (default)\" stringValue=\"basic_status\"/> \
        <String integralValue=\"2\" stringDoc=\"library outputs detailed information\" stringValue=\"detailed_status\"/> \
        <String integralValue=\"3\" stringDoc=\"library outputs very detailed information\" stringValue=\"verbose_detailed_status\"/> \
      </Validator> \
      <Validator defaultParameterName=\"timer_type\" integralValue=\"int\" type=\"StringIntegralValidator(int)\" validatorId=\"2\"> \
        <String integralValue=\"0\" stringDoc=\"No timing data will be collected (the default).\" stringValue=\"no_timers\"/> \
        <String integralValue=\"1\" stringDoc=\"Time an algorithm (or other entity) as a whole.\" stringValue=\"macro_timers\"/> \
        <String integralValue=\"2\" stringDoc=\"Time the substeps of an entity.\" stringValue=\"micro_timers\"/> \
        <String integralValue=\"3\" stringDoc=\"Run both MACRO and MICRO timers.\" stringValue=\"both_timers\"/> \
        <String integralValue=\"4\" stringDoc=\"Run timers added to code for testing, removed later\" stringValue=\"test_timers\"/> \
      </Validator> \
      <Validator defaultParameterName=\"timer_output_stream\" integralValue=\"int\" type=\"StringIntegralValidator(int)\" validatorId=\"4\"> \
        <String integralValue=\"0\" stringValue=\"std::cout\"/> \
        <String integralValue=\"0\" stringValue=\"cout\"/> \
        <String integralValue=\"0\" stringValue=\"stdout\"/> \
        <String integralValue=\"1\" stringValue=\"std::cerr\"/> \
        <String integralValue=\"1\" stringValue=\"cerr\"/> \
        <String integralValue=\"1\" stringValue=\"stderr\"/> \
        <String integralValue=\"2\" stringValue=\"/dev/null\"/> \
        <String integralValue=\"2\" stringValue=\"null\"/> \
      </Validator> \
      <Validator defaultParameterName=\"memory_output_stream\" integralValue=\"int\" type=\"StringIntegralValidator(int)\" validatorId=\"5\"> \
        <String integralValue=\"0\" stringValue=\"std::cout\"/> \
        <String integralValue=\"0\" stringValue=\"cout\"/> \
        <String integralValue=\"0\" stringValue=\"stdout\"/> \
        <String integralValue=\"1\" stringValue=\"std::cerr\"/> \
        <String integralValue=\"1\" stringValue=\"cerr\"/> \
        <String integralValue=\"1\" stringValue=\"stderr\"/> \
        <String integralValue=\"2\" stringValue=\"/dev/null\"/> \
        <String integralValue=\"2\" stringValue=\"null\"/> \
      </Validator> \
      <Validator fileMustExist=\"false\" type=\"FilenameValidator\" validatorId=\"6\"/> \
      <Validator type=\"IntegerRangeListValidator(int)\" unsorted=\"false\" validatorId=\"7\"/> \
      <Validator type=\"IntegerRangeListValidator(int)\" unsorted=\"false\" validatorId=\"9\"/> \
      <Validator type=\"IntegerRangeListValidator(int)\" unsorted=\"true\" validatorId=\"8\"/> \
      <Validator type=\"StringValidator\" validatorId=\"10\"> \
        <String value=\"speed\"/> \
        <String value=\"balance\"/> \
        <String value=\"quality\"/> \
      </Validator> \
      <Validator allowDouble=\"true\" allowInt=\"true\" allowString=\"true\" prefferedType=\"double\" type=\"anynumberValidator\" validatorId=\"12\"/> \
      <Validator type=\"StringValidator\" validatorId=\"11\"> \
        <String value=\"memory\"/> \
        <String value=\"balance\"/> \
        <String value=\"speed\"/> \
      </Validator> \
      <Validator type=\"StringValidator\" validatorId=\"13\"> \
        <String value=\"rcm\"/> \
        <String value=\"minimum_degree\"/> \
        <String value=\"natural\"/> \
        <String value=\"random\"/> \
      </Validator> \
      <Validator defaultParameterName=\"compute_metrics\" integralValue=\"int\" type=\"StringIntegralValidator(int)\" validatorId=\"15\"> \
        <String integralValue=\"1\" stringValue=\"true\"/> \
        <String integralValue=\"1\" stringValue=\"yes\"/> \
        <String integralValue=\"1\" stringValue=\"1\"/> \
        <String integralValue=\"1\" stringValue=\"on\"/> \
        <String integralValue=\"0\" stringValue=\"false\"/> \
        <String integralValue=\"0\" stringValue=\"no\"/> \
        <String integralValue=\"0\" stringValue=\"0\"/> \
        <String integralValue=\"0\" stringValue=\"off\"/> \
      </Validator> \
      <Validator type=\"StringValidator\" validatorId=\"14\"> \
        <String value=\"amd\"/> \
        <String value=\"package2\"/> \
        <String value=\"package3\"/> \
      </Validator> \
      <Validator defaultParameterName=\"randomize_input\" integralValue=\"int\" type=\"StringIntegralValidator(int)\" validatorId=\"17\"> \
        <String integralValue=\"1\" stringValue=\"true\"/> \
        <String integralValue=\"1\" stringValue=\"yes\"/> \
        <String integralValue=\"1\" stringValue=\"1\"/> \
        <String integralValue=\"1\" stringValue=\"on\"/> \
        <String integralValue=\"0\" stringValue=\"false\"/> \
        <String integralValue=\"0\" stringValue=\"no\"/> \
        <String integralValue=\"0\" stringValue=\"0\"/> \
        <String integralValue=\"0\" stringValue=\"off\"/> \
      </Validator> \
      <Validator type=\"IntegerRangeListValidator(int)\" unsorted=\"true\" validatorId=\"16\"/> \
      <Validator allowDouble=\"true\" allowInt=\"true\" allowString=\"true\" prefferedType=\"double\" type=\"anynumberValidator\" validatorId=\"19\"/> \
      <Validator allowDouble=\"true\" allowInt=\"true\" allowString=\"true\" prefferedType=\"double\" type=\"anynumberValidator\" validatorId=\"20\"/> \
      <Validator type=\"StringValidator\" validatorId=\"18\"> \
        <String value=\"balance_object_count\"/> \
        <String value=\"balance_object_weight\"/> \
        <String value=\"multicriteria_minimize_total_weight\"/> \
        <String value=\"multicriteria_minimize_maximum_weight\"/> \
        <String value=\"multicriteria_balance_total_maximum\"/> \
        <String value=\"minimize_cut_edge_count\"/> \
        <String value=\"minimize_cut_edge_weight\"/> \
        <String value=\"minimize_neighboring_parts\"/> \
        <String value=\"minimize_boundary_vertices\"/> \
      </Validator> \
      <Validator allowDouble=\"true\" allowInt=\"true\" allowString=\"true\" prefferedType=\"double\" type=\"anynumberValidator\" validatorId=\"21\"/> \
      <Validator type=\"StringValidator\" validatorId=\"22\"> \
        <String value=\"partition\"/> \
        <String value=\"repartition\"/> \
        <String value=\"maximize_overlap\"/> \
      </Validator> \
      <Validator type=\"StringValidator\" validatorId=\"23\"> \
        <String value=\"matrix_rows\"/> \
        <String value=\"matrix_columns\"/> \
        <String value=\"matrix_nonzeros\"/> \
        <String value=\"mesh_elements\"/> \
        <String value=\"mesh_nodes\"/> \
        <String value=\"graph_edges\"/> \
        <String value=\"graph_vertices\"/> \
        <String value=\"coordinates\"/> \
        <String value=\"identifiers\"/> \
      </Validator> \
      <Validator type=\"StringValidator\" validatorId=\"24\"> \
        <String value=\"hypergraph\"/> \
        <String value=\"graph\"/> \
        <String value=\"geometry\"/> \
        <String value=\"ids\"/> \
      </Validator> \
      <Validator type=\"StringValidator\" validatorId=\"25\"> \
        <String value=\"rcb\"/> \
        <String value=\"multijagged\"/> \
        <String value=\"rib\"/> \
        <String value=\"hsfc\"/> \
        <String value=\"patoh\"/> \
        <String value=\"phg\"/> \
        <String value=\"metis\"/> \
        <String value=\"parmetis\"/> \
        <String value=\"scotch\"/> \
        <String value=\"ptscotch\"/> \
        <String value=\"block\"/> \
        <String value=\"cyclic\"/> \
        <String value=\"random\"/> \
      </Validator> \
      <Validator defaultParameterName=\"rectilinear_blocks\" integralValue=\"int\" type=\"StringIntegralValidator(int)\" validatorId=\"26\"> \
        <String integralValue=\"1\" stringValue=\"true\"/> \
        <String integralValue=\"1\" stringValue=\"yes\"/> \
        <String integralValue=\"1\" stringValue=\"1\"/> \
        <String integralValue=\"1\" stringValue=\"on\"/> \
        <String integralValue=\"0\" stringValue=\"false\"/> \
        <String integralValue=\"0\" stringValue=\"no\"/> \
        <String integralValue=\"0\" stringValue=\"0\"/> \
        <String integralValue=\"0\" stringValue=\"off\"/> \
      </Validator> \
      <Validator defaultParameterName=\"average_cuts\" integralValue=\"int\" type=\"StringIntegralValidator(int)\" validatorId=\"27\"> \
        <String integralValue=\"1\" stringValue=\"true\"/> \
        <String integralValue=\"1\" stringValue=\"yes\"/> \
        <String integralValue=\"1\" stringValue=\"1\"/> \
        <String integralValue=\"1\" stringValue=\"on\"/> \
        <String integralValue=\"0\" stringValue=\"false\"/> \
        <String integralValue=\"0\" stringValue=\"no\"/> \
        <String integralValue=\"0\" stringValue=\"0\"/> \
        <String integralValue=\"0\" stringValue=\"off\"/> \
      </Validator> \
      <Validator max=\"250\" min=\"1\" precision=\"0\" step=\"1\" type=\"EnhancedNumberValidator(int)\" validatorId=\"28\"/> \
      <Validator type=\"StringValidator\" validatorId=\"29\"> \
        <String value=\"no\"/> \
        <String value=\"transpose\"/> \
        <String value=\"bipartite\"/> \
      </Validator> \
      <Validator defaultParameterName=\"subset_graph\" integralValue=\"int\" type=\"StringIntegralValidator(int)\" validatorId=\"30\"> \
        <String integralValue=\"1\" stringValue=\"true\"/> \
        <String integralValue=\"1\" stringValue=\"yes\"/> \
        <String integralValue=\"1\" stringValue=\"1\"/> \
        <String integralValue=\"1\" stringValue=\"on\"/> \
        <String integralValue=\"0\" stringValue=\"false\"/> \
        <String integralValue=\"0\" stringValue=\"no\"/> \
        <String integralValue=\"0\" stringValue=\"0\"/> \
        <String integralValue=\"0\" stringValue=\"off\"/> \
      </Validator> \
      <Validator defaultParameterName=\"parallel_part_calculation_count\" allowDouble=\"false\" allowInt=\"true\" allowString=\"true\" prefferedType=\"int\" type=\"anynumberValidator\" validatorId=\"33\"/> \
      <Validator defaultParameterName=\"migration_all_to_all_type\" allowDouble=\"false\" allowInt=\"true\" allowString=\"true\" prefferedType=\"int\" type=\"anynumberValidator\" validatorId=\"37\"/> \
      <Validator defaultParameterName=\"migration_check_option\" allowDouble=\"false\" allowInt=\"true\" allowString=\"true\" prefferedType=\"int\" type=\"anynumberValidator\" validatorId=\"38\"/> \
      <Validator defaultParameterName=\"migration_processor_assignment_type\" allowDouble=\"false\" allowInt=\"true\" allowString=\"true\" prefferedType=\"int\" type=\"anynumberValidator\" validatorId=\"39\"/> \
      <Validator defaultParameterName=\"remap_parts\" integralValue=\"int\" type=\"StringIntegralValidator(int)\" validatorId=\"40\"> \
        <String integralValue=\"1\" stringValue=\"true\"/> \
        <String integralValue=\"1\" stringValue=\"yes\"/> \
        <String integralValue=\"1\" stringValue=\"1\"/> \
        <String integralValue=\"1\" stringValue=\"on\"/> \
        <String integralValue=\"0\" stringValue=\"false\"/> \
        <String integralValue=\"0\" stringValue=\"no\"/> \
        <String integralValue=\"0\" stringValue=\"0\"/> \
        <String integralValue=\"0\" stringValue=\"off\"/> \
      </Validator> \
    </Validators> \
  </ParameterList>"

#endif  //ZOLTAN2_PARAMETER_DEFINITION_HEADER
