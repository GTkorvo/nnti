#ifndef stk_adapt_IElementBasedAdapterPredicate_hpp
#define stk_adapt_IElementBasedAdapterPredicate_hpp

#include <stk_adapt/IAdapter.hpp>
#include <stk_percept/PerceptMesh.hpp>

#include <functional>

namespace stk {
  namespace adapt {

    /** Signatures for predicate objects that can be used to select entities (elements, edges, faces,...) for
     *   refinement or unrefinement.  The class must supply an operator() that takes an entity and decides
     *   on whether it should be refined, or unrefined, or ignored.
     *
     * The Selector pattern as shown below is useful for selecting only entities that belong to particular
     *   mesh Parts, for example, or any other definition of Selector.
     *
     * We follow the unary_function pattern to enable the structs to be used in STL algorithms that know about
     *   unary_functions.
     *
     * The following are a couple of examples of what refine and unrefine predicates might look like.
     *
     * Following these examples we show the prototype for the operations that are performed on these predicates.
     */

    // Example
    struct IElementBasedAdapterPredicate : public std::unary_function<const stk::mesh::Entity , int> {
      PerceptMesh& m_eMesh;
      stk::mesh::Selector * m_eb_selector;
      stk::mesh::FieldBase *m_field;
      double m_tolerance;
    protected:
      IElementBasedAdapterPredicate(PerceptMesh& eMesh, stk::mesh::Selector * selector=0, stk::mesh::FieldBase *field=0, double tolerance=0.0) :
        m_eMesh(eMesh), m_eb_selector(selector), m_field(field), m_tolerance(tolerance) {}
    };

    // Can be instantiated by the user, or used to define your own
    struct ElementRefinePredicate : public IElementBasedAdapterPredicate {
      //PerceptMesh& m_eMesh;

      ElementRefinePredicate(PerceptMesh& eMesh, stk::mesh::Selector* selector=0, stk::mesh::FieldBase *field=0, double tolerance=0.0) :
        IElementBasedAdapterPredicate(eMesh, selector, field, tolerance) {}

      /// Return DO_REFINE, DO_UNREFINE, DO_NOTHING
      int operator()(const stk::mesh::Entity entity) {
        double *fdata = 0;
        if (m_field)
          fdata = eMesh.field_data( *static_cast<const ScalarFieldType *>(m_field) , entity );
        bool selected = (m_eb_selector==0 || (*m_eb_selector)(m_eMesh.bucket(entity)));
        bool ref_field_criterion = (fdata  && fdata[0] > 0);
        bool unref_field_criterion = (fdata && fdata[0] < 0);
        int mark = 0;
        if (selected && ref_field_criterion) mark |= DO_REFINE;
        if (selected && unref_field_criterion) mark |= DO_UNREFINE;
        return mark;
      }
    };

    /** Example of how the above is used (pseudo-code):
     *
     * class PredicateAdapter : public Refiner {
     *   void visit_for_refine_and_unrefin(std::unary_function<stk::mesh::Entity , bool>& user_predicate_refine)
     *     {
     *        foreach (Entity entity in list)
     *          {
     *             if (user_predicate_refine(entity) & DO_REFINE) { do_refine_operation(entity); }
     *          }
     *        foreach (Entity entity in list)
     *          {
     *             if (user_predicate_unrefine(entity) & DO_UNREFINE) { do_unrefine_operation(entity); }
     *          }
     *     }
     * };
     */

    //get_selected_entities(Selector, rank, EntityVec vec);


    //filter_iterator( Predicate&, vec.begin());
  }
}
#endif
