#ifndef stk_adapt_UniformRefinerPattern_Pyramid13_Pyramid13_10_sierra_hpp
#define stk_adapt_UniformRefinerPattern_Pyramid13_Pyramid13_10_sierra_hpp

#include <stk_adapt/sierra_element/RefinementTopology.hpp>
#include <stk_adapt/sierra_element/StdMeshObjTopologies.hpp>

#include "UniformRefinerPattern_Quad8_Quad8_4_sierra.hpp"
#include "UniformRefinerPattern_Tri6_Tri6_4_sierra.hpp"

#include <stk_percept/PerceptBoostArray.hpp>

namespace stk {
  namespace adapt {

#define DEBUG_Pyramid13_Pyramid13_10 0

    // Some explanation: Pyramid refinement pattern creates a heterogeneous mesh, so we create two
    //   sub-patterns to deal with the two different resulting topologies (pyramid and tet).  A third
    //   (parent) pattern is created to refer to the two sub-patterns, akin to URP_Heterogeneous_3D.

    //================================================================================================================================================================
    template <>
    class UniformRefinerPattern<shards::Pyramid<13>, shards::Pyramid<13>, 6, SierraPort > : public URP<shards::Pyramid<13>, shards::Pyramid<13>  >
    {


    public:

      UniformRefinerPattern(percept::PerceptMesh& eMesh, BlockNamesType block_names = BlockNamesType()) :  URP<shards::Pyramid<13>, shards::Pyramid<13>  >(eMesh)
      {
        m_primaryEntityRank = stk::mesh::MetaData::ELEMENT_RANK;
        Elem::StdMeshObjTopologies::bootstrap();
        m_do_strip_hashes=false;
      }

      ~UniformRefinerPattern()
      {
      }

      void setSubPatterns( std::vector<UniformRefinerPatternBase *>& bp, percept::PerceptMesh& eMesh )
      {
        EXCEPTWATCH;
        bp = std::vector<UniformRefinerPatternBase *>(1u, 0);

        bp[0] = this;
      }

      virtual void doBreak() {}
      void fillNeededEntities(std::vector<NeededEntityType>& needed_entities)
      {
        needed_entities.resize(3);
        needed_entities[0] = NeededEntityType(m_eMesh.edge_rank(), 3u);
        needed_entities[1] = NeededEntityType(m_eMesh.face_rank(), 9u);   // cheating here - we re-use the full quadratic face
        needed_entities[2] = NeededEntityType(stk::mesh::MetaData::ELEMENT_RANK, 4u);
      }

      virtual unsigned getNumNewElemPerElem() { return 6; }


      void
      createNewElements(percept::PerceptMesh& eMesh, NodeRegistry& nodeRegistry,
                        stk::mesh::Entity element,  NewSubEntityNodesType& new_sub_entity_nodes, vector<stk::mesh::Entity>::iterator& element_pool,
                        stk::mesh::FieldBase *proc_rank_field=0)
      {
        genericRefine_createNewElements(eMesh, nodeRegistry,
                                        element, new_sub_entity_nodes, element_pool,
                                        proc_rank_field);
      }

    };

    //================================================================================================================================================================
    template <>
    class UniformRefinerPattern<shards::Pyramid<13>, shards::Tetrahedron<10>, 4, SierraPort > : public URP<shards::Pyramid<13>, shards::Tetrahedron<10>  >
    {

    public:

      UniformRefinerPattern(percept::PerceptMesh& eMesh, BlockNamesType block_names = BlockNamesType()) :  URP<shards::Pyramid<13>, shards::Tetrahedron<10>  >(eMesh)
      {
        m_primaryEntityRank = stk::mesh::MetaData::ELEMENT_RANK;

        Elem::StdMeshObjTopologies::bootstrap();
        m_do_strip_hashes=false;

      }

      ~UniformRefinerPattern()
      {
      }

      void setSubPatterns( std::vector<UniformRefinerPatternBase *>& bp, percept::PerceptMesh& eMesh )
      {
        EXCEPTWATCH;
        bp = std::vector<UniformRefinerPatternBase *>(1u, 0);

        bp[0] = this;
      }

      virtual void doBreak() {}
      void fillNeededEntities(std::vector<NeededEntityType>& needed_entities)
      {
        needed_entities.resize(3);
        needed_entities[0] = NeededEntityType(m_eMesh.edge_rank(), 3u);
        needed_entities[1] = NeededEntityType(m_eMesh.face_rank(), 9u);   // cheating here - we re-use the full quadratic face
        needed_entities[2] = NeededEntityType(stk::mesh::MetaData::ELEMENT_RANK, 4u);
      }

      virtual unsigned getNumNewElemPerElem() { return 4; }

      void
      createNewElements(percept::PerceptMesh& eMesh, NodeRegistry& nodeRegistry,
                        stk::mesh::Entity element,  NewSubEntityNodesType& new_sub_entity_nodes, vector<stk::mesh::Entity>::iterator& element_pool,
                        stk::mesh::FieldBase *proc_rank_field=0)
      {
        genericRefine_createNewElements(eMesh, nodeRegistry,
                                        element, new_sub_entity_nodes, element_pool,
                                        proc_rank_field);
      }

    };

    //================================================================================================================================================================
    template <>
    class UniformRefinerPattern<shards::Pyramid<13>, shards::Pyramid<13>, 10, SierraPort > : public URP<shards::Pyramid<13>, shards::Pyramid<13>  >
                               //class UniformRefinerPattern<shards::Pyramid<13>, shards::Pyramid<13>, 10, SierraPort > : public UniformRefinerPatternBase
    {
      UniformRefinerPattern<shards::Quadrilateral<8>, shards::Quadrilateral<8>, 4, SierraPort > * m_face_breaker;
      UniformRefinerPattern<shards::Triangle<6>, shards::Triangle<6>, 4, SierraPort > * m_face_breaker_tri;

      std::vector<UniformRefinerPatternBase *> m_bp;
    protected:

      percept::PerceptMesh& m_eMesh;

    public:

      static void printParts(UniformRefinerPatternBase *bp)
      {
        for (unsigned ii=0; ii < bp->getFromParts().size(); ii++)
          {
            std::cout << "tmp Pyramid13_Pyramid13_10 ii, fromParts= " << ii << " " << bp->getFromParts()[ii]->name() << std::endl;
          }
        for (unsigned ii=0; ii < bp->getToParts().size(); ii++)
          {
            std::cout << "tmp Pyramid13_Pyramid13_10 ii, toParts= " << ii << " " << bp->getToParts()[ii]->name() << std::endl;
          }
      }

      //UniformRefinerPattern(percept::PerceptMesh& eMesh, BlockNamesType block_names = BlockNamesType()) :  URP<shards::Pyramid<13>, shards::Pyramid<13>  >(eMesh), m_eMesh(eMesh)
      UniformRefinerPattern(percept::PerceptMesh& eMesh, BlockNamesType block_names = BlockNamesType()) : URP<shards::Pyramid<13>, shards::Pyramid<13>  >(eMesh), m_eMesh(eMesh)
      {
        m_primaryEntityRank = stk::mesh::MetaData::ELEMENT_RANK;
        m_do_strip_hashes=false;

        Elem::StdMeshObjTopologies::bootstrap();

        // list all types of known break patterns to be used here
        m_bp.resize(0);

        m_bp.push_back(new  UniformRefinerPattern<shards::Pyramid<13>,       shards::Pyramid<13>,      6, SierraPort > (eMesh, block_names));
        m_bp.push_back(new  UniformRefinerPattern<shards::Pyramid<13>,       shards::Tetrahedron<10>,  4, SierraPort > (eMesh, block_names));

        bool sameTopology = false;
        //setNeededParts(eMesh, block_names, sameTopology);
        m_bp[0]->setNeededParts(eMesh, block_names, sameTopology); // force a new part for pyramids
        if (DEBUG_Pyramid13_Pyramid13_10)
          {
            std::cout << "tmp Pyramid13_Pyramid13_10 printParts m_bp[0]= " ; printParts(m_bp[0]);
          }
        m_bp[1]->setNeededParts(eMesh, block_names, sameTopology);
        if (DEBUG_Pyramid13_Pyramid13_10)
          {
            std::cout << "tmp Pyramid13_Pyramid13_10 printParts m_bp[1]= " ; printParts(m_bp[1]);
          }

        for (int ibp=0; ibp < 2; ibp++)
          {
            stk::mesh::PartVector& fromParts = m_bp[ibp]->getFromParts();
            for (unsigned ii=0; ii < fromParts.size(); ii++)
              {
                if (std::find(getFromParts().begin(), getFromParts().end(), fromParts[ii]) == getFromParts().end())
                  {
                    getFromParts().push_back(fromParts[ii]);
                  }
              }
            stk::mesh::PartVector& toParts = m_bp[ibp]->getToParts();
            for (unsigned ii=0; ii < toParts.size(); ii++)
              {
                if (std::find(getToParts().begin(), getToParts().end(), toParts[ii]) == getToParts().end())
                  {
                    getToParts().push_back(toParts[ii]);
                  }
              }
          }
        if (DEBUG_Pyramid13_Pyramid13_10)
          {
            std::cout << "tmp Pyramid13_Pyramid13_10 printParts this= " ;
            printParts(this);
          }

        m_face_breaker =  new UniformRefinerPattern<shards::Quadrilateral<8>, shards::Quadrilateral<8>, 4, SierraPort > (eMesh, block_names) ;
        m_face_breaker_tri = new UniformRefinerPattern<shards::Triangle<6>, shards::Triangle<6>, 4, SierraPort > (eMesh, block_names);
      }

      ~UniformRefinerPattern()
      {
        if (m_face_breaker) delete m_face_breaker;
        if (m_face_breaker_tri) delete m_face_breaker_tri;
        for (unsigned ibp=0; ibp < m_bp.size(); ibp++)
          {
            if (m_bp[ibp]) delete m_bp[ibp];
          }
      }

      void setSubPatterns( std::vector<UniformRefinerPatternBase *>& bp, percept::PerceptMesh& eMesh )
      {
        EXCEPTWATCH;

        bp.resize(3u);
        bp[0] = this;
        bp[1] = m_face_breaker;
        bp[2] = m_face_breaker_tri;
      }

      virtual void doBreak()
      {
        throw std::runtime_error("shouldn't call URP_Pyramid13_Pyramid13::doBreak()");
      }
      virtual unsigned getFromTypeKey()
      {
        //throw std::runtime_error("shouldn't call URP_Pyramid13_Pyramid13::getFromTypeKey()");
        return shards::Pyramid<13>::key;
      }

      virtual std::string getFromTopoPartName() {
        shards::CellTopology cell_topo(getFromTopology());
        return cell_topo.getName();
      }
      virtual std::string getToTopoPartName() {
        shards::CellTopology cell_topo(getToTopology());
        return cell_topo.getName();
      }

      virtual const CellTopologyData * getFromTopology() { return shards::getCellTopologyData< shards::Pyramid<13> >(); }
      virtual const CellTopologyData * getToTopology() { return shards::getCellTopologyData< shards::Pyramid<13> >(); }

      //       virtual const CellTopologyData *  getFromTopology()
      //       {
      //         throw std::runtime_error("shouldn't call URP_Pyramid13_Pyramid13::getFromTopology()");
      //       }

      //       virtual const CellTopologyData *  getToTopology() {
      //         throw std::runtime_error("shouldn't call URP_Pyramid13_Pyramid13::getToTopology()");
      //       }

      void fillNeededEntities(std::vector<NeededEntityType>& needed_entities)
      {
        //throw std::runtime_error("shouldn't call URP_Pyramid13_Pyramid13::fillNeededEntities()");
        needed_entities.resize(3);
        needed_entities[0] = NeededEntityType(m_eMesh.edge_rank(), 3u);
        needed_entities[1] = NeededEntityType(m_eMesh.face_rank(), 9u);   // cheating here - we re-use the full quadratic face
        needed_entities[2] = NeededEntityType(stk::mesh::MetaData::ELEMENT_RANK, 4u);
      }

      virtual unsigned getNumNewElemPerElem()
      {
        //throw std::runtime_error("shouldn't call URP_Pyramid13_Pyramid13::getNumNewElemPerElem()");
        return 10;
      }

      void
      createNewElements(percept::PerceptMesh& eMesh, NodeRegistry& nodeRegistry,
                        stk::mesh::Entity element,  NewSubEntityNodesType& new_sub_entity_nodes, vector<stk::mesh::Entity>::iterator& element_pool,
                        stk::mesh::FieldBase *proc_rank_field=0)
      {
        //throw std::runtime_error("shouldn't call URP_Pyramid13_Pyramid13::createNewElements()");

        // pyramids
        m_bp[0]->createNewElements(eMesh, nodeRegistry, element, new_sub_entity_nodes, element_pool, proc_rank_field);
        // tets
        m_bp[1]->createNewElements(eMesh, nodeRegistry, element, new_sub_entity_nodes, element_pool, proc_rank_field);
      }

    };

  }
}
#endif
