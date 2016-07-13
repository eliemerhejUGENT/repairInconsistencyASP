//
//  main.cpp
//  repairInconsistentASP
//
//  Created by Elie Merhej on 01/03/15.
//  Copyright (c) 2015 Elie Merhej. All rights reserved.
//

#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <cmath>

enum EDGE_TYPE
{
    ACTIVATES = 0,
    INHIBITS,
};

enum TABLE_TYPE
{
    ACTIVE = 0,
    INACTIVE,
};

enum
{
    NOT_CORRUPTED = 0,
    CORRUPTED = 1,
};

struct Edge
{
    Edge(int Type, unsigned int From, unsigned int To)
    {
        type = Type;
        from = From;
        to = To;
    }
    
    bool operator==(const Edge &rhs)
    {
        return ((type == rhs.type) && (from == rhs.from) && (to == rhs.to));
    }
    
    int type;
    unsigned int from;
    unsigned int to;
};

struct TableElement
{
    TableElement(int Type, unsigned int Gene, unsigned int Time)
    {
        type = Type;
        gene = Gene;
        time = Time;
    }
    
    int type;
    unsigned int gene;
    unsigned int time;
};

struct GeneNetwork
{
    GeneNetwork()
    {
        timeSteps = 0;
        geneNum = 0;
        
        kDegree = 0.0f;
    }
    
    GeneNetwork(GeneNetwork& rhs)
    {
        name = rhs.name;
        
        timeSteps = rhs.timeSteps;
        for(size_t i = 0; i < rhs.table.size(); ++i)
            table.push_back(rhs.table[i]);
        
        geneNum = rhs.geneNum;
        
        for(size_t i = 0; i < rhs.edges.size(); ++i)
            edges.push_back(rhs.edges[i]);
        
        for(size_t i = 0; i < rhs.addedEdges.size(); ++i)
            addedEdges.push_back(rhs.addedEdges[i]);
        
        kDegree = rhs.kDegree;
        edgesNodesRatio = rhs.edgesNodesRatio;
        diameter = rhs.diameter;
        
    }
    
    void LearnProperties()
    {
        if(edges.empty())
        {
            std::cout << "\nERROR: Network not loaded, can't learn properties...\n\n";
            return;
        }
        
        CalculateKDegree();
        CalculateEdgesNodesRatio();
    }
    
    void CalculateKDegree()
    {
        unsigned int kDegrees = 0;
        
        for(size_t i = 1; i <= geneNum; ++i)
        {
            unsigned int geneKDegree = 0;
            
            for(size_t j = 0; j < edges.size(); ++j)
            {
                if(edges[j].from == i)
                    ++geneKDegree;
                
                if(edges[j].to == i)
                    ++geneKDegree;
            }
            
            kDegrees += geneKDegree;
        }
        
        kDegree = (float)kDegrees / (float)geneNum;
    }
    
    void CalculateEdgesNodesRatio()
    {
        edgesNodesRatio = (float)edges.size() / (float)geneNum;
    }
    
    void PrintProperties()
    {
        if(edges.empty())
        {
            std::cout << "\nERROR: Network not loaded, can't print properties...\n\n";
            return;
        }
        
        std::cout << "\n\nNetwork: " << name << "\n";
        std::cout << "\nNb. of nodes = " << geneNum << "\n";
        std::cout << "\nNb. of edges = " << edges.size() << "\n";
        std::cout << "\nEdges-Nodes Ratio = " << edgesNodesRatio << "\n";
        std::cout << "\nAverage K degree of all nodes = " << kDegree << "\n";
        std::cout << "\nDiameter = " << diameter << "\n\n\n";
    }
    
    std::string name;
    
    unsigned int timeSteps;
    std::vector<TableElement> table;
    
    unsigned int geneNum;
    std::vector<Edge> edges;
    std::vector<Edge> addedEdges; // edges added to corrupt the network
    
    float kDegree; // average of number of total edges connected to a node
    float edgesNodesRatio; // ratio of edges per node
    unsigned int diameter; // diameter of the network (largest value of smallest distances between every pair of nodes)
};




struct Repair
{
    Repair(unsigned int r0 = 0, unsigned int r1 = 0, unsigned int r2 = 0, unsigned int r3 = 0, unsigned int r4 = 0, unsigned int r5 = 0, unsigned int r6 = 0)
    {
        ruleViolations[0] = r0;
        ruleViolations[1] = r1;
        ruleViolations[2] = r2;
        ruleViolations[3] = r3;
        ruleViolations[4] = r4;
        ruleViolations[5] = r5;
        ruleViolations[6] = r6;
        
        for(size_t i = 0; i < 7; ++i)
            zScores[i] = 0.0f;
        
        totalZScore = 0.0f;
    }
    
    unsigned int ruleViolations[7];
    float zScores[7];
    
    float totalZScore;
};







GeneNetwork budding;
GeneNetwork fission;
GeneNetwork elegans;
GeneNetwork mammalian;
GeneNetwork arabidopsis;
GeneNetwork thcell;







void LoadNetworks(int Status)
{
    //***********************************************************************************************************************************
    //***********************************************************************************************************************************
    //***********************************************************************************************************************************
    budding.name = "budding";
    budding.timeSteps = 13;
    budding.geneNum = 11;
    budding.diameter = 3; // learned from separate ASP program
    
    if(Status == NOT_CORRUPTED)
    {
        // Edges of original network
        budding.edges.push_back(Edge(ACTIVATES,1,3));
        budding.edges.push_back(Edge(ACTIVATES,1,2));
        budding.edges.push_back(Edge(INHIBITS,1,1));
        budding.edges.push_back(Edge(ACTIVATES,3,4));
        budding.edges.push_back(Edge(ACTIVATES,2,8));
        budding.edges.push_back(Edge(INHIBITS,4,9));
        budding.edges.push_back(Edge(INHIBITS,4,5));
        budding.edges.push_back(Edge(INHIBITS,4,4));
        budding.edges.push_back(Edge(ACTIVATES,8,11));
        budding.edges.push_back(Edge(ACTIVATES,8,10));
        budding.edges.push_back(Edge(INHIBITS,8,9));
        budding.edges.push_back(Edge(INHIBITS,8,5));
        budding.edges.push_back(Edge(INHIBITS,9,10));
        budding.edges.push_back(Edge(INHIBITS,9,8));
        budding.edges.push_back(Edge(INHIBITS,5,10));
        budding.edges.push_back(Edge(ACTIVATES,11,10));
        budding.edges.push_back(Edge(ACTIVATES,11,7));
        budding.edges.push_back(Edge(ACTIVATES,11,6));
        budding.edges.push_back(Edge(INHIBITS,11,11));
        budding.edges.push_back(Edge(ACTIVATES,10,11));
        budding.edges.push_back(Edge(ACTIVATES,10,7));
        budding.edges.push_back(Edge(INHIBITS,10,3));
        budding.edges.push_back(Edge(INHIBITS,10,2));
        budding.edges.push_back(Edge(INHIBITS,10,9));
        budding.edges.push_back(Edge(INHIBITS,10,5));
        budding.edges.push_back(Edge(INHIBITS,10,6));
        budding.edges.push_back(Edge(ACTIVATES,7,5));
        budding.edges.push_back(Edge(ACTIVATES,7,9));
        budding.edges.push_back(Edge(ACTIVATES,7,6));
        budding.edges.push_back(Edge(INHIBITS,7,10));
        budding.edges.push_back(Edge(INHIBITS,7,8));
        budding.edges.push_back(Edge(INHIBITS,7,7));
        budding.edges.push_back(Edge(ACTIVATES,6,9));
        budding.edges.push_back(Edge(INHIBITS,6,6));
    }
    else if(Status == CORRUPTED)
    {
        // Edges of original network
        budding.edges.push_back(Edge(ACTIVATES,1,3));
        budding.edges.push_back(Edge(ACTIVATES,1,2));
        budding.edges.push_back(Edge(INHIBITS,1,1));
        budding.edges.push_back(Edge(ACTIVATES,3,4));
        //budding.edges.push_back(Edge(ACTIVATES,2,8));
        budding.edges.push_back(Edge(INHIBITS,4,9));
        budding.edges.push_back(Edge(INHIBITS,4,5));
        budding.edges.push_back(Edge(INHIBITS,4,4));
        budding.edges.push_back(Edge(ACTIVATES,8,11));
        budding.edges.push_back(Edge(ACTIVATES,8,10));
        budding.edges.push_back(Edge(INHIBITS,8,9));
        budding.edges.push_back(Edge(INHIBITS,8,5));
        budding.edges.push_back(Edge(INHIBITS,9,10));
        //budding.edges.push_back(Edge(INHIBITS,9,8));
        budding.edges.push_back(Edge(INHIBITS,5,10));
        budding.edges.push_back(Edge(ACTIVATES,11,10));
        budding.edges.push_back(Edge(ACTIVATES,11,7));
        budding.edges.push_back(Edge(ACTIVATES,11,6));
        budding.edges.push_back(Edge(INHIBITS,11,11));
        //budding.edges.push_back(Edge(ACTIVATES,10,11));
        budding.edges.push_back(Edge(ACTIVATES,10,7));
        //budding.edges.push_back(Edge(INHIBITS,10,3));
        //budding.edges.push_back(Edge(INHIBITS,10,2));
        //budding.edges.push_back(Edge(INHIBITS,10,9));
        budding.edges.push_back(Edge(INHIBITS,10,5));
        budding.edges.push_back(Edge(INHIBITS,10,6));
        budding.edges.push_back(Edge(ACTIVATES,7,5));
        budding.edges.push_back(Edge(ACTIVATES,7,9));
        budding.edges.push_back(Edge(ACTIVATES,7,6));
        budding.edges.push_back(Edge(INHIBITS,7,10));
        //budding.edges.push_back(Edge(INHIBITS,7,8));
        budding.edges.push_back(Edge(INHIBITS,7,7));
        budding.edges.push_back(Edge(ACTIVATES,6,9));
        budding.edges.push_back(Edge(INHIBITS,6,6));
        
        
        budding.addedEdges.push_back(Edge(ACTIVATES,1,5));
        budding.addedEdges.push_back(Edge(ACTIVATES,1,8));
        budding.addedEdges.push_back(Edge(ACTIVATES,1,9));
        budding.addedEdges.push_back(Edge(ACTIVATES,1,11));
        budding.addedEdges.push_back(Edge(ACTIVATES,2,1));
        budding.addedEdges.push_back(Edge(ACTIVATES,2,4));
        budding.addedEdges.push_back(Edge(ACTIVATES,3,8));
        budding.addedEdges.push_back(Edge(ACTIVATES,7,3));
        budding.addedEdges.push_back(Edge(ACTIVATES,7,11));
        budding.addedEdges.push_back(Edge(ACTIVATES,9,3));
        budding.addedEdges.push_back(Edge(INHIBITS,3,1));
        budding.addedEdges.push_back(Edge(INHIBITS,5,8));
        budding.addedEdges.push_back(Edge(INHIBITS,6,1));
        budding.addedEdges.push_back(Edge(INHIBITS,6,2));
        budding.addedEdges.push_back(Edge(INHIBITS,6,3));
        budding.addedEdges.push_back(Edge(INHIBITS,6,4));
        budding.addedEdges.push_back(Edge(INHIBITS,6,5));
        budding.addedEdges.push_back(Edge(INHIBITS,6,7));
        budding.addedEdges.push_back(Edge(INHIBITS,6,8));
        budding.addedEdges.push_back(Edge(INHIBITS,6,10));
        budding.addedEdges.push_back(Edge(INHIBITS,7,1));
        budding.addedEdges.push_back(Edge(INHIBITS,7,2));
        budding.addedEdges.push_back(Edge(INHIBITS,8,8));
        budding.addedEdges.push_back(Edge(INHIBITS,9,11));
        budding.addedEdges.push_back(Edge(INHIBITS,11,2));
        budding.addedEdges.push_back(Edge(INHIBITS,11,3));
        budding.addedEdges.push_back(Edge(INHIBITS,11,5));
        budding.addedEdges.push_back(Edge(INHIBITS,11,9));
    }
    
    // Timeseries table
    budding.table.push_back(TableElement(ACTIVE,1,1));
    budding.table.push_back(TableElement(INACTIVE,2,1));
    budding.table.push_back(TableElement(INACTIVE,3,1));
    budding.table.push_back(TableElement(INACTIVE,4,1));
    budding.table.push_back(TableElement(ACTIVE,5,1));
    budding.table.push_back(TableElement(INACTIVE,6,1));
    budding.table.push_back(TableElement(INACTIVE,7,1));
    budding.table.push_back(TableElement(INACTIVE,8,1));
    budding.table.push_back(TableElement(ACTIVE,9,1));
    budding.table.push_back(TableElement(INACTIVE,10,1));
    budding.table.push_back(TableElement(INACTIVE,11,1));
    
    budding.table.push_back(TableElement(INACTIVE,1,2));
    budding.table.push_back(TableElement(ACTIVE,2,2));
    budding.table.push_back(TableElement(ACTIVE,3,2));
    budding.table.push_back(TableElement(INACTIVE,4,2));
    budding.table.push_back(TableElement(ACTIVE,5,2));
    budding.table.push_back(TableElement(INACTIVE,6,2));
    budding.table.push_back(TableElement(INACTIVE,7,2));
    budding.table.push_back(TableElement(INACTIVE,8,2));
    budding.table.push_back(TableElement(ACTIVE,9,2));
    budding.table.push_back(TableElement(INACTIVE,10,2));
    budding.table.push_back(TableElement(INACTIVE,11,2));
    
    budding.table.push_back(TableElement(INACTIVE,1,3));
    budding.table.push_back(TableElement(ACTIVE,2,3));
    budding.table.push_back(TableElement(ACTIVE,3,3));
    budding.table.push_back(TableElement(ACTIVE,4,3));
    budding.table.push_back(TableElement(ACTIVE,5,3));
    budding.table.push_back(TableElement(INACTIVE,6,3));
    budding.table.push_back(TableElement(INACTIVE,7,3));
    budding.table.push_back(TableElement(INACTIVE,8,3));
    budding.table.push_back(TableElement(ACTIVE,9,3));
    budding.table.push_back(TableElement(INACTIVE,10,3));
    budding.table.push_back(TableElement(INACTIVE,11,3));
    
    budding.table.push_back(TableElement(INACTIVE,1,4));
    budding.table.push_back(TableElement(ACTIVE,2,4));
    budding.table.push_back(TableElement(ACTIVE,3,4));
    budding.table.push_back(TableElement(ACTIVE,4,4));
    budding.table.push_back(TableElement(INACTIVE,5,4));
    budding.table.push_back(TableElement(INACTIVE,6,4));
    budding.table.push_back(TableElement(INACTIVE,7,4));
    budding.table.push_back(TableElement(INACTIVE,8,4));
    budding.table.push_back(TableElement(INACTIVE,9,4));
    budding.table.push_back(TableElement(INACTIVE,10,4));
    budding.table.push_back(TableElement(INACTIVE,11,4));
    
    budding.table.push_back(TableElement(INACTIVE,1,5));
    budding.table.push_back(TableElement(ACTIVE,2,5));
    budding.table.push_back(TableElement(ACTIVE,3,5));
    budding.table.push_back(TableElement(ACTIVE,4,5));
    budding.table.push_back(TableElement(INACTIVE,5,5));
    budding.table.push_back(TableElement(INACTIVE,6,5));
    budding.table.push_back(TableElement(INACTIVE,7,5));
    budding.table.push_back(TableElement(ACTIVE,8,5));
    budding.table.push_back(TableElement(INACTIVE,9,5));
    budding.table.push_back(TableElement(INACTIVE,10,5));
    budding.table.push_back(TableElement(INACTIVE,11,5));
    
    budding.table.push_back(TableElement(INACTIVE,1,6));
    budding.table.push_back(TableElement(ACTIVE,2,6));
    budding.table.push_back(TableElement(ACTIVE,3,6));
    budding.table.push_back(TableElement(ACTIVE,4,6));
    budding.table.push_back(TableElement(INACTIVE,5,6));
    budding.table.push_back(TableElement(INACTIVE,6,6));
    budding.table.push_back(TableElement(INACTIVE,7,6));
    budding.table.push_back(TableElement(ACTIVE,8,6));
    budding.table.push_back(TableElement(INACTIVE,9,6));
    budding.table.push_back(TableElement(ACTIVE,10,6));
    budding.table.push_back(TableElement(ACTIVE,11,6));
    
    budding.table.push_back(TableElement(INACTIVE,1,7));
    budding.table.push_back(TableElement(INACTIVE,2,7));
    budding.table.push_back(TableElement(INACTIVE,3,7));
    budding.table.push_back(TableElement(ACTIVE,4,7));
    budding.table.push_back(TableElement(INACTIVE,5,7));
    budding.table.push_back(TableElement(INACTIVE,6,7));
    budding.table.push_back(TableElement(ACTIVE,7,7));
    budding.table.push_back(TableElement(ACTIVE,8,7));
    budding.table.push_back(TableElement(INACTIVE,9,7));
    budding.table.push_back(TableElement(ACTIVE,10,7));
    budding.table.push_back(TableElement(ACTIVE,11,7));
    
    budding.table.push_back(TableElement(INACTIVE,1,8));
    budding.table.push_back(TableElement(INACTIVE,2,8));
    budding.table.push_back(TableElement(INACTIVE,3,8));
    budding.table.push_back(TableElement(INACTIVE,4,8));
    budding.table.push_back(TableElement(INACTIVE,5,8));
    budding.table.push_back(TableElement(INACTIVE,6,8));
    budding.table.push_back(TableElement(ACTIVE,7,8));
    budding.table.push_back(TableElement(INACTIVE,8,8));
    budding.table.push_back(TableElement(INACTIVE,9,8));
    budding.table.push_back(TableElement(ACTIVE,10,8));
    budding.table.push_back(TableElement(ACTIVE,11,8));
    
    budding.table.push_back(TableElement(INACTIVE,1,9));
    budding.table.push_back(TableElement(INACTIVE,2,9));
    budding.table.push_back(TableElement(INACTIVE,3,9));
    budding.table.push_back(TableElement(INACTIVE,4,9));
    budding.table.push_back(TableElement(INACTIVE,5,9));
    budding.table.push_back(TableElement(INACTIVE,6,9));
    budding.table.push_back(TableElement(ACTIVE,7,9));
    budding.table.push_back(TableElement(INACTIVE,8,9));
    budding.table.push_back(TableElement(INACTIVE,9,9));
    budding.table.push_back(TableElement(ACTIVE,10,9));
    budding.table.push_back(TableElement(ACTIVE,11,9));
    
    budding.table.push_back(TableElement(INACTIVE,1,10));
    budding.table.push_back(TableElement(INACTIVE,2,10));
    budding.table.push_back(TableElement(INACTIVE,3,10));
    budding.table.push_back(TableElement(INACTIVE,4,10));
    budding.table.push_back(TableElement(INACTIVE,5,10));
    budding.table.push_back(TableElement(INACTIVE,6,10));
    budding.table.push_back(TableElement(ACTIVE,7,10));
    budding.table.push_back(TableElement(INACTIVE,8,10));
    budding.table.push_back(TableElement(INACTIVE,9,10));
    budding.table.push_back(TableElement(ACTIVE,10,10));
    budding.table.push_back(TableElement(ACTIVE,11,10));
    
    budding.table.push_back(TableElement(INACTIVE,1,11));
    budding.table.push_back(TableElement(INACTIVE,2,11));
    budding.table.push_back(TableElement(INACTIVE,3,11));
    budding.table.push_back(TableElement(INACTIVE,4,11));
    budding.table.push_back(TableElement(INACTIVE,5,11));
    budding.table.push_back(TableElement(INACTIVE,6,11));
    budding.table.push_back(TableElement(ACTIVE,7,11));
    budding.table.push_back(TableElement(INACTIVE,8,11));
    budding.table.push_back(TableElement(INACTIVE,9,11));
    budding.table.push_back(TableElement(ACTIVE,10,11));
    budding.table.push_back(TableElement(ACTIVE,11,11));
    
    budding.table.push_back(TableElement(INACTIVE,1,12));
    budding.table.push_back(TableElement(INACTIVE,2,12));
    budding.table.push_back(TableElement(INACTIVE,3,12));
    budding.table.push_back(TableElement(INACTIVE,4,12));
    budding.table.push_back(TableElement(INACTIVE,5,12));
    budding.table.push_back(TableElement(INACTIVE,6,12));
    budding.table.push_back(TableElement(ACTIVE,7,12));
    budding.table.push_back(TableElement(INACTIVE,8,12));
    budding.table.push_back(TableElement(INACTIVE,9,12));
    budding.table.push_back(TableElement(ACTIVE,10,12));
    budding.table.push_back(TableElement(ACTIVE,11,12));
    
    budding.table.push_back(TableElement(INACTIVE,1,13));
    budding.table.push_back(TableElement(INACTIVE,2,13));
    budding.table.push_back(TableElement(INACTIVE,3,13));
    budding.table.push_back(TableElement(INACTIVE,4,13));
    budding.table.push_back(TableElement(INACTIVE,5,13));
    budding.table.push_back(TableElement(INACTIVE,6,13));
    budding.table.push_back(TableElement(ACTIVE,7,13));
    budding.table.push_back(TableElement(INACTIVE,8,13));
    budding.table.push_back(TableElement(INACTIVE,9,13));
    budding.table.push_back(TableElement(ACTIVE,10,13));
    budding.table.push_back(TableElement(ACTIVE,11,13));
    
    //safety check
    size_t buddingTableSize = budding.geneNum * budding.timeSteps;
    size_t buddingTableElementsNum = budding.table.size();
    
    if(buddingTableSize != buddingTableElementsNum)
    {
        std::cout << "\nERROR: missing entries in Budding timeseries table...\n\n";
        return;
    }
    
    //***********************************************************************************************************************************
    //***********************************************************************************************************************************
    //***********************************************************************************************************************************
    fission.name = "fission";
    fission.timeSteps = 13;
    fission.geneNum = 9;
    fission.diameter = 3; // learned from separate ASP program
    
    if(Status == NOT_CORRUPTED)
    {
        // Edges of original network
        fission.edges.push_back(Edge(INHIBITS,1,3));
        fission.edges.push_back(Edge(INHIBITS,1,2));
        fission.edges.push_back(Edge(INHIBITS,1,1));
        fission.edges.push_back(Edge(INHIBITS,2,6));
        fission.edges.push_back(Edge(INHIBITS,2,9));
        fission.edges.push_back(Edge(INHIBITS,3,6));
        fission.edges.push_back(Edge(INHIBITS,3,9));
        fission.edges.push_back(Edge(ACTIVATES,4,3));
        fission.edges.push_back(Edge(ACTIVATES,4,2));
        fission.edges.push_back(Edge(ACTIVATES,4,8));
        fission.edges.push_back(Edge(INHIBITS,4,5));
        fission.edges.push_back(Edge(INHIBITS,4,4));
        fission.edges.push_back(Edge(ACTIVATES,5,9));
        fission.edges.push_back(Edge(ACTIVATES,6,5));
        fission.edges.push_back(Edge(INHIBITS,6,3));
        fission.edges.push_back(Edge(INHIBITS,6,2));
        fission.edges.push_back(Edge(INHIBITS,6,8));
        fission.edges.push_back(Edge(ACTIVATES,7,4));
        fission.edges.push_back(Edge(INHIBITS,7,6));
        fission.edges.push_back(Edge(INHIBITS,7,9));
        fission.edges.push_back(Edge(INHIBITS,7,7));
        fission.edges.push_back(Edge(INHIBITS,8,9));
        fission.edges.push_back(Edge(ACTIVATES,9,7));
        fission.edges.push_back(Edge(INHIBITS,9,3));
        fission.edges.push_back(Edge(INHIBITS,9,2));
    }
    else if(Status == CORRUPTED)
    {
        fission.edges.push_back(Edge(INHIBITS,1,3));
        fission.edges.push_back(Edge(INHIBITS,1,2));
        fission.edges.push_back(Edge(INHIBITS,1,1));
        //fission.edges.push_back(Edge(INHIBITS,2,6));
        //fission.edges.push_back(Edge(INHIBITS,2,9));
        //fission.edges.push_back(Edge(INHIBITS,3,6));
        //fission.edges.push_back(Edge(INHIBITS,3,9));
        fission.edges.push_back(Edge(ACTIVATES,4,3));
        fission.edges.push_back(Edge(ACTIVATES,4,2));
        fission.edges.push_back(Edge(ACTIVATES,4,8));
        fission.edges.push_back(Edge(INHIBITS,4,5));
        fission.edges.push_back(Edge(INHIBITS,4,4));
        fission.edges.push_back(Edge(ACTIVATES,5,9));
        fission.edges.push_back(Edge(ACTIVATES,6,5));
        fission.edges.push_back(Edge(INHIBITS,6,3));
        fission.edges.push_back(Edge(INHIBITS,6,2));
        fission.edges.push_back(Edge(INHIBITS,6,8));
        fission.edges.push_back(Edge(ACTIVATES,7,4));
        fission.edges.push_back(Edge(INHIBITS,7,6));
        fission.edges.push_back(Edge(INHIBITS,7,9));
        fission.edges.push_back(Edge(INHIBITS,7,7));
        //fission.edges.push_back(Edge(INHIBITS,8,9));
        fission.edges.push_back(Edge(ACTIVATES,9,7));
        fission.edges.push_back(Edge(INHIBITS,9,3));
        fission.edges.push_back(Edge(INHIBITS,9,2));
        
        fission.addedEdges.push_back(Edge(ACTIVATES,2,2));
        fission.addedEdges.push_back(Edge(ACTIVATES,2,3));
        fission.addedEdges.push_back(Edge(ACTIVATES,2,8));
        fission.addedEdges.push_back(Edge(ACTIVATES,3,2));
        fission.addedEdges.push_back(Edge(ACTIVATES,3,3));
        fission.addedEdges.push_back(Edge(ACTIVATES,3,8));
        fission.addedEdges.push_back(Edge(ACTIVATES,8,2));
        fission.addedEdges.push_back(Edge(ACTIVATES,8,3));
        fission.addedEdges.push_back(Edge(ACTIVATES,8,8));
        fission.addedEdges.push_back(Edge(INHIBITS,3,4));
        fission.addedEdges.push_back(Edge(INHIBITS,3,5));
        fission.addedEdges.push_back(Edge(INHIBITS,3,7));
        fission.addedEdges.push_back(Edge(INHIBITS,4,7));
        fission.addedEdges.push_back(Edge(INHIBITS,7,2));
        fission.addedEdges.push_back(Edge(INHIBITS,7,3));
        fission.addedEdges.push_back(Edge(INHIBITS,8,1));
        fission.addedEdges.push_back(Edge(INHIBITS,8,5));
        fission.addedEdges.push_back(Edge(INHIBITS,8,6));
        fission.addedEdges.push_back(Edge(INHIBITS,9,4));
        fission.addedEdges.push_back(Edge(INHIBITS,9,9));
        
    }
    
    // Timeseries table
    fission.table.push_back(TableElement(ACTIVE,1,1));
    fission.table.push_back(TableElement(ACTIVE,4,1));
    fission.table.push_back(TableElement(ACTIVE,5,1));
    fission.table.push_back(TableElement(INACTIVE,2,1));
    fission.table.push_back(TableElement(INACTIVE,3,1));
    fission.table.push_back(TableElement(INACTIVE,6,1));
    fission.table.push_back(TableElement(INACTIVE,7,1));
    fission.table.push_back(TableElement(INACTIVE,8,1));
    fission.table.push_back(TableElement(INACTIVE,9,1));
    fission.table.push_back(TableElement(INACTIVE,2,2));
    fission.table.push_back(TableElement(INACTIVE,7,2));
    fission.table.push_back(TableElement(INACTIVE,6,2));
    fission.table.push_back(TableElement(INACTIVE,3,2));
    fission.table.push_back(TableElement(INACTIVE,2,3));
    fission.table.push_back(TableElement(INACTIVE,9,3));
    fission.table.push_back(TableElement(INACTIVE,6,3));
    fission.table.push_back(TableElement(INACTIVE,3,3));
    fission.table.push_back(TableElement(INACTIVE,1,2));
    fission.table.push_back(TableElement(INACTIVE,5,2));
    fission.table.push_back(TableElement(INACTIVE,4,2));
    fission.table.push_back(TableElement(ACTIVE,4,4));
    fission.table.push_back(TableElement(ACTIVE,8,2));
    fission.table.push_back(TableElement(ACTIVE,9,2));
    fission.table.push_back(TableElement(INACTIVE,2,4));
    fission.table.push_back(TableElement(INACTIVE,9,4));
    fission.table.push_back(TableElement(INACTIVE,7,4));
    fission.table.push_back(TableElement(INACTIVE,6,4));
    fission.table.push_back(TableElement(INACTIVE,3,4));
    fission.table.push_back(TableElement(INACTIVE,1,3));
    fission.table.push_back(TableElement(INACTIVE,5,3));
    fission.table.push_back(TableElement(INACTIVE,4,3));
    fission.table.push_back(TableElement(ACTIVE,8,3));
    fission.table.push_back(TableElement(INACTIVE,9,5));
    fission.table.push_back(TableElement(INACTIVE,7,5));
    fission.table.push_back(TableElement(INACTIVE,6,5));
    fission.table.push_back(TableElement(INACTIVE,1,4));
    fission.table.push_back(TableElement(INACTIVE,5,4));
    fission.table.push_back(TableElement(ACTIVE,8,4));
    fission.table.push_back(TableElement(INACTIVE,9,6));
    fission.table.push_back(TableElement(INACTIVE,7,6));
    fission.table.push_back(TableElement(INACTIVE,6,6));
    fission.table.push_back(TableElement(INACTIVE,1,5));
    fission.table.push_back(TableElement(INACTIVE,5,5));
    fission.table.push_back(TableElement(INACTIVE,4,5));
    fission.table.push_back(TableElement(ACTIVE,2,5));
    fission.table.push_back(TableElement(ACTIVE,8,5));
    fission.table.push_back(TableElement(ACTIVE,3,5));
    fission.table.push_back(TableElement(ACTIVE,7,3));
    fission.table.push_back(TableElement(INACTIVE,9,7));
    fission.table.push_back(TableElement(INACTIVE,7,7));
    fission.table.push_back(TableElement(INACTIVE,6,7));
    fission.table.push_back(TableElement(INACTIVE,1,6));
    fission.table.push_back(TableElement(INACTIVE,5,6));
    fission.table.push_back(TableElement(INACTIVE,4,6));
    fission.table.push_back(TableElement(ACTIVE,2,6));
    fission.table.push_back(TableElement(ACTIVE,8,6));
    fission.table.push_back(TableElement(ACTIVE,3,6));
    fission.table.push_back(TableElement(INACTIVE,9,8));
    fission.table.push_back(TableElement(INACTIVE,7,8));
    fission.table.push_back(TableElement(INACTIVE,6,8));
    fission.table.push_back(TableElement(INACTIVE,1,7));
    fission.table.push_back(TableElement(INACTIVE,5,7));
    fission.table.push_back(TableElement(INACTIVE,4,7));
    fission.table.push_back(TableElement(ACTIVE,2,7));
    fission.table.push_back(TableElement(ACTIVE,8,7));
    fission.table.push_back(TableElement(ACTIVE,3,7));
    fission.table.push_back(TableElement(INACTIVE,9,9));
    fission.table.push_back(TableElement(INACTIVE,7,9));
    fission.table.push_back(TableElement(INACTIVE,6,9));
    fission.table.push_back(TableElement(INACTIVE,1,8));
    fission.table.push_back(TableElement(INACTIVE,5,8));
    fission.table.push_back(TableElement(INACTIVE,4,8));
    fission.table.push_back(TableElement(ACTIVE,2,8));
    fission.table.push_back(TableElement(ACTIVE,8,8));
    fission.table.push_back(TableElement(ACTIVE,3,8));
    fission.table.push_back(TableElement(INACTIVE,9,10));
    fission.table.push_back(TableElement(INACTIVE,7,10));
    fission.table.push_back(TableElement(INACTIVE,6,10));
    fission.table.push_back(TableElement(INACTIVE,1,9));
    fission.table.push_back(TableElement(INACTIVE,5,9));
    fission.table.push_back(TableElement(INACTIVE,4,9));
    fission.table.push_back(TableElement(ACTIVE,2,9));
    fission.table.push_back(TableElement(ACTIVE,8,9));
    fission.table.push_back(TableElement(ACTIVE,3,9));
    fission.table.push_back(TableElement(INACTIVE,9,11));
    fission.table.push_back(TableElement(INACTIVE,7,11));
    fission.table.push_back(TableElement(INACTIVE,6,11));
    fission.table.push_back(TableElement(INACTIVE,1,10));
    fission.table.push_back(TableElement(INACTIVE,5,10));
    fission.table.push_back(TableElement(INACTIVE,4,10));
    fission.table.push_back(TableElement(ACTIVE,2,10));
    fission.table.push_back(TableElement(ACTIVE,8,10));
    fission.table.push_back(TableElement(ACTIVE,3,10));
    fission.table.push_back(TableElement(INACTIVE,9,12));
    fission.table.push_back(TableElement(INACTIVE,7,12));
    fission.table.push_back(TableElement(INACTIVE,6,12));
    fission.table.push_back(TableElement(INACTIVE,1,11));
    fission.table.push_back(TableElement(INACTIVE,5,11));
    fission.table.push_back(TableElement(INACTIVE,4,11));
    fission.table.push_back(TableElement(ACTIVE,2,11));
    fission.table.push_back(TableElement(ACTIVE,8,11));
    fission.table.push_back(TableElement(ACTIVE,3,11));
    fission.table.push_back(TableElement(INACTIVE,9,13));
    fission.table.push_back(TableElement(INACTIVE,7,13));
    fission.table.push_back(TableElement(INACTIVE,6,13));
    fission.table.push_back(TableElement(INACTIVE,1,12));
    fission.table.push_back(TableElement(INACTIVE,5,12));
    fission.table.push_back(TableElement(INACTIVE,4,12));
    fission.table.push_back(TableElement(ACTIVE,2,12));
    fission.table.push_back(TableElement(ACTIVE,8,12));
    fission.table.push_back(TableElement(ACTIVE,3,12));
    fission.table.push_back(TableElement(INACTIVE,1,13));
    fission.table.push_back(TableElement(INACTIVE,5,13));
    fission.table.push_back(TableElement(INACTIVE,4,13));
    fission.table.push_back(TableElement(ACTIVE,2,13));
    fission.table.push_back(TableElement(ACTIVE,8,13));
    fission.table.push_back(TableElement(ACTIVE,3,13));
    
    //safety check
    size_t fissionTableSize = fission.geneNum * fission.timeSteps;
    size_t fissionTableElementsNum = fission.table.size();
    
    if(fissionTableSize != fissionTableElementsNum)
    {
        std::cout << "\nERROR: missing entries in Fission timeseries table...\n\n";
        return;
    }
    
    
    //***********************************************************************************************************************************
    //***********************************************************************************************************************************
    //***********************************************************************************************************************************
    elegans.name = "elegans";
    elegans.timeSteps = 7;
    elegans.geneNum = 8;
    elegans.diameter = 3; // learned from separate ASP program
    
    if(Status == NOT_CORRUPTED)
    {
        // Edges of original network
        elegans.edges.push_back(Edge(ACTIVATES,1,7));
        elegans.edges.push_back(Edge(INHIBITS,1,6));
        elegans.edges.push_back(Edge(INHIBITS,1,3));
        elegans.edges.push_back(Edge(INHIBITS,2,1));
        elegans.edges.push_back(Edge(INHIBITS,2,2));
        elegans.edges.push_back(Edge(INHIBITS,2,5));
        elegans.edges.push_back(Edge(ACTIVATES,3,2));
        elegans.edges.push_back(Edge(INHIBITS,3,5));
        elegans.edges.push_back(Edge(ACTIVATES,3,4));
        elegans.edges.push_back(Edge(ACTIVATES,4,2));
        elegans.edges.push_back(Edge(INHIBITS,4,5));
        elegans.edges.push_back(Edge(INHIBITS,4,4));
        elegans.edges.push_back(Edge(ACTIVATES,5,3));
        elegans.edges.push_back(Edge(ACTIVATES,6,5));
        elegans.edges.push_back(Edge(INHIBITS,6,1));
        elegans.edges.push_back(Edge(ACTIVATES,7,6));
        elegans.edges.push_back(Edge(INHIBITS,7,1));
        elegans.edges.push_back(Edge(INHIBITS,7,8));
        elegans.edges.push_back(Edge(INHIBITS,7,7));
        elegans.edges.push_back(Edge(ACTIVATES,8,5));
        elegans.edges.push_back(Edge(ACTIVATES,8,1));
    }
    else if(Status == CORRUPTED)
    {
        elegans.edges.push_back(Edge(ACTIVATES,1,7));
        elegans.edges.push_back(Edge(INHIBITS,1,6));
        elegans.edges.push_back(Edge(INHIBITS,1,3));
        //elegans.edges.push_back(Edge(INHIBITS,2,1));
        elegans.edges.push_back(Edge(INHIBITS,2,2));
        elegans.edges.push_back(Edge(INHIBITS,2,5));
        elegans.edges.push_back(Edge(ACTIVATES,3,2));
        //elegans.edges.push_back(Edge(INHIBITS,3,5));
        elegans.edges.push_back(Edge(ACTIVATES,3,4));
        elegans.edges.push_back(Edge(ACTIVATES,4,2));
        elegans.edges.push_back(Edge(INHIBITS,4,5));
        //elegans.edges.push_back(Edge(INHIBITS,4,4));
        elegans.edges.push_back(Edge(ACTIVATES,5,3));
        elegans.edges.push_back(Edge(ACTIVATES,6,5));
        //elegans.edges.push_back(Edge(INHIBITS,6,1));
        elegans.edges.push_back(Edge(ACTIVATES,7,6));
        elegans.edges.push_back(Edge(INHIBITS,7,1));
        //elegans.edges.push_back(Edge(INHIBITS,7,8));
        elegans.edges.push_back(Edge(INHIBITS,7,7));
        elegans.edges.push_back(Edge(ACTIVATES,8,5));
        elegans.edges.push_back(Edge(ACTIVATES,8,1));
        
        elegans.addedEdges.push_back(Edge(ACTIVATES,3,7));
        elegans.addedEdges.push_back(Edge(ACTIVATES,5,6));
        elegans.addedEdges.push_back(Edge(ACTIVATES,7,2));
        elegans.addedEdges.push_back(Edge(ACTIVATES,7,3));
        elegans.addedEdges.push_back(Edge(ACTIVATES,7,5));
        elegans.addedEdges.push_back(Edge(INHIBITS,1,8));
        elegans.addedEdges.push_back(Edge(INHIBITS,2,3));
        elegans.addedEdges.push_back(Edge(INHIBITS,2,8));
        elegans.addedEdges.push_back(Edge(INHIBITS,3,1));
        elegans.addedEdges.push_back(Edge(INHIBITS,3,8));
        elegans.addedEdges.push_back(Edge(INHIBITS,4,8));
        elegans.addedEdges.push_back(Edge(INHIBITS,5,1));
        elegans.addedEdges.push_back(Edge(INHIBITS,5,5));
        elegans.addedEdges.push_back(Edge(INHIBITS,5,7));
        elegans.addedEdges.push_back(Edge(INHIBITS,5,8));
        elegans.addedEdges.push_back(Edge(INHIBITS,6,8));
        elegans.addedEdges.push_back(Edge(INHIBITS,7,4));
        
    }
    
    
    // Timeseries table
    elegans.table.push_back(TableElement(ACTIVE,1,1));
    elegans.table.push_back(TableElement(ACTIVE,2,1));
    elegans.table.push_back(TableElement(ACTIVE,3,1));
    elegans.table.push_back(TableElement(ACTIVE,4,1));
    elegans.table.push_back(TableElement(ACTIVE,6,1));
    elegans.table.push_back(TableElement(INACTIVE,5,1));
    elegans.table.push_back(TableElement(INACTIVE,7,1));
    elegans.table.push_back(TableElement(INACTIVE,8,1));
    elegans.table.push_back(TableElement(ACTIVE,2,2));
    elegans.table.push_back(TableElement(ACTIVE,4,2));
    elegans.table.push_back(TableElement(INACTIVE,5,2));
    elegans.table.push_back(TableElement(INACTIVE,8,2));
    elegans.table.push_back(TableElement(ACTIVE,2,3));
    elegans.table.push_back(TableElement(ACTIVE,6,3));
    elegans.table.push_back(TableElement(INACTIVE,5,3));
    elegans.table.push_back(TableElement(INACTIVE,7,3));
    elegans.table.push_back(TableElement(INACTIVE,8,3));
    elegans.table.push_back(TableElement(INACTIVE,1,2));
    elegans.table.push_back(TableElement(INACTIVE,3,2));
    elegans.table.push_back(TableElement(INACTIVE,6,2));
    elegans.table.push_back(TableElement(ACTIVE,6,4));
    elegans.table.push_back(TableElement(ACTIVE,7,2));
    elegans.table.push_back(TableElement(INACTIVE,5,4));
    elegans.table.push_back(TableElement(INACTIVE,7,4));
    elegans.table.push_back(TableElement(INACTIVE,8,4));
    elegans.table.push_back(TableElement(INACTIVE,1,3));
    elegans.table.push_back(TableElement(INACTIVE,3,3));
    elegans.table.push_back(TableElement(INACTIVE,4,3));
    elegans.table.push_back(TableElement(ACTIVE,6,5));
    elegans.table.push_back(TableElement(INACTIVE,7,5));
    elegans.table.push_back(TableElement(INACTIVE,8,5));
    elegans.table.push_back(TableElement(INACTIVE,1,4));
    elegans.table.push_back(TableElement(INACTIVE,2,4));
    elegans.table.push_back(TableElement(INACTIVE,3,4));
    elegans.table.push_back(TableElement(INACTIVE,4,4));
    elegans.table.push_back(TableElement(ACTIVE,3,6));
    elegans.table.push_back(TableElement(ACTIVE,6,6));
    elegans.table.push_back(TableElement(INACTIVE,7,6));
    elegans.table.push_back(TableElement(INACTIVE,8,6));
    elegans.table.push_back(TableElement(INACTIVE,1,5));
    elegans.table.push_back(TableElement(INACTIVE,2,5));
    elegans.table.push_back(TableElement(INACTIVE,3,5));
    elegans.table.push_back(TableElement(INACTIVE,4,5));
    elegans.table.push_back(TableElement(ACTIVE,2,7));
    elegans.table.push_back(TableElement(ACTIVE,3,7));
    elegans.table.push_back(TableElement(ACTIVE,4,7));
    elegans.table.push_back(TableElement(ACTIVE,5,5));
    elegans.table.push_back(TableElement(ACTIVE,6,7));
    elegans.table.push_back(TableElement(INACTIVE,7,7));
    elegans.table.push_back(TableElement(INACTIVE,8,7));
    elegans.table.push_back(TableElement(INACTIVE,1,6));
    elegans.table.push_back(TableElement(INACTIVE,2,6));
    elegans.table.push_back(TableElement(INACTIVE,4,6));
    elegans.table.push_back(TableElement(ACTIVE,5,6));
    elegans.table.push_back(TableElement(INACTIVE,1,7));
    elegans.table.push_back(TableElement(ACTIVE,5,7));
    
    
    
    //safety check
    size_t elegansTableSize = elegans.geneNum * elegans.timeSteps;
    size_t elegansTableElementsNum = elegans.table.size();
    
    if(elegansTableSize != elegansTableElementsNum)
    {
        std::cout << "\nERROR: missing entries in Elegans timeseries table...\n\n";
        return;
    }
    
    
    //***********************************************************************************************************************************
    //***********************************************************************************************************************************
    //***********************************************************************************************************************************
    mammalian.name = "mammalian";
    mammalian.timeSteps = 5;
    mammalian.geneNum = 10;
    mammalian.diameter = 3; // learned from separate ASP program
    
    if(Status == NOT_CORRUPTED)
    {
        // Edges of original network
        mammalian.edges.push_back(Edge(ACTIVATES,1,1));
        mammalian.edges.push_back(Edge(INHIBITS,1,2));
        mammalian.edges.push_back(Edge(INHIBITS,1,3));
        mammalian.edges.push_back(Edge(INHIBITS,2,4));
        mammalian.edges.push_back(Edge(INHIBITS,2,5));
        mammalian.edges.push_back(Edge(INHIBITS,2,6));
        mammalian.edges.push_back(Edge(ACTIVATES,3,3));
        mammalian.edges.push_back(Edge(ACTIVATES,3,4));
        mammalian.edges.push_back(Edge(ACTIVATES,3,5));
        mammalian.edges.push_back(Edge(ACTIVATES,3,9));
        mammalian.edges.push_back(Edge(ACTIVATES,3,2));
        mammalian.edges.push_back(Edge(INHIBITS,4,3));
        mammalian.edges.push_back(Edge(INHIBITS,4,4));
        mammalian.edges.push_back(Edge(INHIBITS,4,2));
        mammalian.edges.push_back(Edge(ACTIVATES,5,4));
        mammalian.edges.push_back(Edge(ACTIVATES,5,6));
        mammalian.edges.push_back(Edge(ACTIVATES,6,6));
        mammalian.edges.push_back(Edge(INHIBITS,6,3));
        mammalian.edges.push_back(Edge(INHIBITS,6,9));
        mammalian.edges.push_back(Edge(ACTIVATES,6,8));
        mammalian.edges.push_back(Edge(INHIBITS,6,2));
        mammalian.edges.push_back(Edge(INHIBITS,6,5));
        mammalian.edges.push_back(Edge(INHIBITS,6,4));
        mammalian.edges.push_back(Edge(ACTIVATES,7,9));
        mammalian.edges.push_back(Edge(INHIBITS,7,10));
        mammalian.edges.push_back(Edge(ACTIVATES,7,8));
        mammalian.edges.push_back(Edge(INHIBITS,7,6));
        mammalian.edges.push_back(Edge(INHIBITS,8,6));
        mammalian.edges.push_back(Edge(ACTIVATES,8,8));
        mammalian.edges.push_back(Edge(INHIBITS,9,6));
        mammalian.edges.push_back(Edge(INHIBITS,9,7));
        mammalian.edges.push_back(Edge(INHIBITS,9,8));
        mammalian.edges.push_back(Edge(INHIBITS,9,10));
        mammalian.edges.push_back(Edge(ACTIVATES,10,7));
        mammalian.edges.push_back(Edge(INHIBITS,10,9));
        mammalian.edges.push_back(Edge(INHIBITS,10,3));
        mammalian.edges.push_back(Edge(INHIBITS,10,2));
        mammalian.edges.push_back(Edge(INHIBITS,10,5));
        mammalian.edges.push_back(Edge(ACTIVATES,10,8));
    }
    else if(Status == CORRUPTED)
    {
        mammalian.edges.push_back(Edge(ACTIVATES,1,1));
        mammalian.edges.push_back(Edge(INHIBITS,1,2));
        mammalian.edges.push_back(Edge(INHIBITS,1,3));
        mammalian.edges.push_back(Edge(INHIBITS,2,4));
        mammalian.edges.push_back(Edge(INHIBITS,2,5));
        mammalian.edges.push_back(Edge(INHIBITS,2,6));
        mammalian.edges.push_back(Edge(ACTIVATES,3,3));
        //mammalian.edges.push_back(Edge(ACTIVATES,3,4));
        //mammalian.edges.push_back(Edge(ACTIVATES,3,5));
        //mammalian.edges.push_back(Edge(ACTIVATES,3,9));
        mammalian.edges.push_back(Edge(ACTIVATES,3,2));
        mammalian.edges.push_back(Edge(INHIBITS,4,3));
        mammalian.edges.push_back(Edge(INHIBITS,4,4));
        mammalian.edges.push_back(Edge(INHIBITS,4,2));
        mammalian.edges.push_back(Edge(ACTIVATES,5,4));
        mammalian.edges.push_back(Edge(ACTIVATES,5,6));
        mammalian.edges.push_back(Edge(ACTIVATES,6,6));
        mammalian.edges.push_back(Edge(INHIBITS,6,3));
        mammalian.edges.push_back(Edge(INHIBITS,6,9));
        mammalian.edges.push_back(Edge(ACTIVATES,6,8));
        mammalian.edges.push_back(Edge(INHIBITS,6,2));
        mammalian.edges.push_back(Edge(INHIBITS,6,5));
        //mammalian.edges.push_back(Edge(INHIBITS,6,4));
        mammalian.edges.push_back(Edge(ACTIVATES,7,9));
        //mammalian.edges.push_back(Edge(INHIBITS,7,10));
        mammalian.edges.push_back(Edge(ACTIVATES,7,8));
        mammalian.edges.push_back(Edge(INHIBITS,7,6));
        mammalian.edges.push_back(Edge(INHIBITS,8,6));
        //mammalian.edges.push_back(Edge(ACTIVATES,8,8));
        mammalian.edges.push_back(Edge(INHIBITS,9,6));
        mammalian.edges.push_back(Edge(INHIBITS,9,7));
        mammalian.edges.push_back(Edge(INHIBITS,9,8));
        //mammalian.edges.push_back(Edge(INHIBITS,9,10));
        mammalian.edges.push_back(Edge(ACTIVATES,10,7));
        mammalian.edges.push_back(Edge(INHIBITS,10,9));
        //mammalian.edges.push_back(Edge(INHIBITS,10,3));
        mammalian.edges.push_back(Edge(INHIBITS,10,2));
        mammalian.edges.push_back(Edge(INHIBITS,10,5));
        mammalian.edges.push_back(Edge(ACTIVATES,10,8));
        
        mammalian.addedEdges.push_back(Edge(ACTIVATES,1,5));
        mammalian.addedEdges.push_back(Edge(ACTIVATES,1,6));
        mammalian.addedEdges.push_back(Edge(ACTIVATES,1,7));
        mammalian.addedEdges.push_back(Edge(ACTIVATES,1,8));
        mammalian.addedEdges.push_back(Edge(ACTIVATES,1,9));
        mammalian.addedEdges.push_back(Edge(ACTIVATES,1,10));
        mammalian.addedEdges.push_back(Edge(ACTIVATES,2,1));
        mammalian.addedEdges.push_back(Edge(ACTIVATES,3,8));
        mammalian.addedEdges.push_back(Edge(ACTIVATES,4,5));
        mammalian.addedEdges.push_back(Edge(ACTIVATES,4,8));
        mammalian.addedEdges.push_back(Edge(ACTIVATES,4,9));
        mammalian.addedEdges.push_back(Edge(ACTIVATES,5,8));
        mammalian.addedEdges.push_back(Edge(ACTIVATES,7,3));
        mammalian.addedEdges.push_back(Edge(ACTIVATES,7,5));
        mammalian.addedEdges.push_back(Edge(ACTIVATES,8,3));
        mammalian.addedEdges.push_back(Edge(ACTIVATES,8,4));
        mammalian.addedEdges.push_back(Edge(ACTIVATES,9,3));
        mammalian.addedEdges.push_back(Edge(ACTIVATES,10,1));
        mammalian.addedEdges.push_back(Edge(ACTIVATES,10,4));
        mammalian.addedEdges.push_back(Edge(INHIBITS,2,7));
        mammalian.addedEdges.push_back(Edge(INHIBITS,2,8));
        mammalian.addedEdges.push_back(Edge(INHIBITS,2,9));
        mammalian.addedEdges.push_back(Edge(INHIBITS,2,10));
        mammalian.addedEdges.push_back(Edge(INHIBITS,3,1));
        mammalian.addedEdges.push_back(Edge(INHIBITS,4,6));
        mammalian.addedEdges.push_back(Edge(INHIBITS,5,7));
        mammalian.addedEdges.push_back(Edge(INHIBITS,6,1));
        mammalian.addedEdges.push_back(Edge(INHIBITS,6,7));
        mammalian.addedEdges.push_back(Edge(INHIBITS,6,10));
        mammalian.addedEdges.push_back(Edge(INHIBITS,7,2));
        mammalian.addedEdges.push_back(Edge(INHIBITS,8,10));
        mammalian.addedEdges.push_back(Edge(INHIBITS,9,2));
        
        
    }
    
    
    // Timeseries table
    mammalian.table.push_back(TableElement(INACTIVE,1,1));
    mammalian.table.push_back(TableElement(INACTIVE,2,1));
    mammalian.table.push_back(TableElement(INACTIVE,4,1));
    mammalian.table.push_back(TableElement(INACTIVE,5,1));
    mammalian.table.push_back(TableElement(INACTIVE,6,1));
    mammalian.table.push_back(TableElement(INACTIVE,7,1));
    mammalian.table.push_back(TableElement(INACTIVE,8,1));
    mammalian.table.push_back(TableElement(INACTIVE,9,1));
    mammalian.table.push_back(TableElement(ACTIVE,3,1));
    mammalian.table.push_back(TableElement(ACTIVE,10,1));
    mammalian.table.push_back(TableElement(ACTIVE,3,2));
    mammalian.table.push_back(TableElement(ACTIVE,10,2));
    mammalian.table.push_back(TableElement(INACTIVE,1,2));
    mammalian.table.push_back(TableElement(INACTIVE,2,2));
    mammalian.table.push_back(TableElement(INACTIVE,5,2));
    mammalian.table.push_back(TableElement(INACTIVE,6,2));
    mammalian.table.push_back(TableElement(INACTIVE,9,2));
    mammalian.table.push_back(TableElement(ACTIVE,3,3));
    mammalian.table.push_back(TableElement(INACTIVE,1,3));
    mammalian.table.push_back(TableElement(INACTIVE,2,3));
    mammalian.table.push_back(TableElement(INACTIVE,5,3));
    mammalian.table.push_back(TableElement(INACTIVE,6,3));
    mammalian.table.push_back(TableElement(INACTIVE,9,3));
    mammalian.table.push_back(TableElement(ACTIVE,3,4));
    mammalian.table.push_back(TableElement(ACTIVE,4,2));
    mammalian.table.push_back(TableElement(ACTIVE,7,2));
    mammalian.table.push_back(TableElement(ACTIVE,8,2));
    mammalian.table.push_back(TableElement(INACTIVE,1,4));
    mammalian.table.push_back(TableElement(INACTIVE,2,4));
    mammalian.table.push_back(TableElement(INACTIVE,6,4));
    mammalian.table.push_back(TableElement(ACTIVE,3,5));
    mammalian.table.push_back(TableElement(ACTIVE,4,3));
    mammalian.table.push_back(TableElement(ACTIVE,7,3));
    mammalian.table.push_back(TableElement(ACTIVE,8,3));
    mammalian.table.push_back(TableElement(INACTIVE,1,5));
    mammalian.table.push_back(TableElement(INACTIVE,2,5));
    mammalian.table.push_back(TableElement(INACTIVE,6,5));
    mammalian.table.push_back(TableElement(INACTIVE,7,5));
    mammalian.table.push_back(TableElement(ACTIVE,4,4));
    mammalian.table.push_back(TableElement(ACTIVE,5,4));
    mammalian.table.push_back(TableElement(ACTIVE,7,4));
    mammalian.table.push_back(TableElement(ACTIVE,8,4));
    mammalian.table.push_back(TableElement(ACTIVE,9,4));
    mammalian.table.push_back(TableElement(INACTIVE,10,3));
    mammalian.table.push_back(TableElement(ACTIVE,4,5));
    mammalian.table.push_back(TableElement(ACTIVE,5,5));
    mammalian.table.push_back(TableElement(ACTIVE,8,5));
    mammalian.table.push_back(TableElement(ACTIVE,9,5));
    mammalian.table.push_back(TableElement(INACTIVE,10,4));
    mammalian.table.push_back(TableElement(INACTIVE,10,5));
    
    
    
    
    //safety check
    size_t mammalianTableSize = mammalian.geneNum * mammalian.timeSteps;
    size_t mammalianTableElementsNum = mammalian.table.size();
    
    if(mammalianTableSize != mammalianTableElementsNum)
    {
        std::cout << "\nERROR: missing entries in Mammalian timeseries table...\n\n";
        return;
    }
    
    
    //***********************************************************************************************************************************
    //***********************************************************************************************************************************
    //***********************************************************************************************************************************
    arabidopsis.name = "arabidopsis";
    arabidopsis.timeSteps = 5;
    arabidopsis.geneNum = 10;
    arabidopsis.diameter = 4; // learned from separate ASP program
    
    if(Status == NOT_CORRUPTED)
    {
        // Edges of original network
        arabidopsis.edges.push_back(Edge(ACTIVATES,1,3));
        arabidopsis.edges.push_back(Edge(INHIBITS,1,2));
        arabidopsis.edges.push_back(Edge(INHIBITS,1,6));
        arabidopsis.edges.push_back(Edge(INHIBITS,2,4));
        arabidopsis.edges.push_back(Edge(ACTIVATES,2,6));
        arabidopsis.edges.push_back(Edge(INHIBITS,3,4));
        arabidopsis.edges.push_back(Edge(INHIBITS,3,6));
        arabidopsis.edges.push_back(Edge(INHIBITS,4,2));
        arabidopsis.edges.push_back(Edge(INHIBITS,5,4));
        arabidopsis.edges.push_back(Edge(ACTIVATES,6,4));
        arabidopsis.edges.push_back(Edge(ACTIVATES,6,2));
        arabidopsis.edges.push_back(Edge(INHIBITS,6,3));
        arabidopsis.edges.push_back(Edge(ACTIVATES,6,9));
        arabidopsis.edges.push_back(Edge(ACTIVATES,6,7));
        arabidopsis.edges.push_back(Edge(ACTIVATES,7,7));
        arabidopsis.edges.push_back(Edge(ACTIVATES,7,9));
        arabidopsis.edges.push_back(Edge(INHIBITS,8,9));
        arabidopsis.edges.push_back(Edge(INHIBITS,8,7));
        arabidopsis.edges.push_back(Edge(ACTIVATES,9,9));
        arabidopsis.edges.push_back(Edge(ACTIVATES,9,7));
        arabidopsis.edges.push_back(Edge(ACTIVATES,10,7));
        arabidopsis.edges.push_back(Edge(ACTIVATES,10,9));
    }
    else if(Status == CORRUPTED)
    {
        //arabidopsis.edges.push_back(Edge(ACTIVATES,1,3));
        arabidopsis.edges.push_back(Edge(INHIBITS,1,2));
        //arabidopsis.edges.push_back(Edge(INHIBITS,1,6));
        arabidopsis.edges.push_back(Edge(INHIBITS,2,4));
        arabidopsis.edges.push_back(Edge(ACTIVATES,2,6));
        arabidopsis.edges.push_back(Edge(INHIBITS,3,4));
        arabidopsis.edges.push_back(Edge(INHIBITS,3,6));
        arabidopsis.edges.push_back(Edge(INHIBITS,4,2));
        arabidopsis.edges.push_back(Edge(INHIBITS,5,4));
        arabidopsis.edges.push_back(Edge(ACTIVATES,6,4));
        arabidopsis.edges.push_back(Edge(ACTIVATES,6,2));
        arabidopsis.edges.push_back(Edge(INHIBITS,6,3));
        //arabidopsis.edges.push_back(Edge(ACTIVATES,6,9));
        //arabidopsis.edges.push_back(Edge(ACTIVATES,6,7));
        arabidopsis.edges.push_back(Edge(ACTIVATES,7,7));
        arabidopsis.edges.push_back(Edge(ACTIVATES,7,9));
        arabidopsis.edges.push_back(Edge(INHIBITS,8,9));
        arabidopsis.edges.push_back(Edge(INHIBITS,8,7));
        arabidopsis.edges.push_back(Edge(ACTIVATES,9,9));
        //arabidopsis.edges.push_back(Edge(ACTIVATES,9,7));
        arabidopsis.edges.push_back(Edge(ACTIVATES,10,7));
        arabidopsis.edges.push_back(Edge(ACTIVATES,10,9));
        
        arabidopsis.addedEdges.push_back(Edge(ACTIVATES,1,7));
        arabidopsis.addedEdges.push_back(Edge(ACTIVATES,1,9));
        arabidopsis.addedEdges.push_back(Edge(ACTIVATES,4,3));
        arabidopsis.addedEdges.push_back(Edge(ACTIVATES,7,1));
        arabidopsis.addedEdges.push_back(Edge(ACTIVATES,10,6));
        arabidopsis.addedEdges.push_back(Edge(ACTIVATES,10,8));
        arabidopsis.addedEdges.push_back(Edge(ACTIVATES,10,10));
        arabidopsis.addedEdges.push_back(Edge(INHIBITS,2,8));
        arabidopsis.addedEdges.push_back(Edge(INHIBITS,2,9));
        arabidopsis.addedEdges.push_back(Edge(INHIBITS,2,10));
        arabidopsis.addedEdges.push_back(Edge(INHIBITS,5,3));
        arabidopsis.addedEdges.push_back(Edge(INHIBITS,5,5));
        arabidopsis.addedEdges.push_back(Edge(INHIBITS,5,6));
        arabidopsis.addedEdges.push_back(Edge(INHIBITS,5,7));
        arabidopsis.addedEdges.push_back(Edge(INHIBITS,5,8));
        arabidopsis.addedEdges.push_back(Edge(INHIBITS,5,10));
        arabidopsis.addedEdges.push_back(Edge(INHIBITS,6,6));
        arabidopsis.addedEdges.push_back(Edge(INHIBITS,7,6));
        
        
    }
    
    
    // Timeseries table
    arabidopsis.table.push_back(TableElement(ACTIVE,1,1));
    arabidopsis.table.push_back(TableElement(ACTIVE,6,1));
    arabidopsis.table.push_back(TableElement(INACTIVE,2,1));
    arabidopsis.table.push_back(TableElement(INACTIVE,3,1));
    arabidopsis.table.push_back(TableElement(INACTIVE,4,1));
    arabidopsis.table.push_back(TableElement(INACTIVE,5,1));
    arabidopsis.table.push_back(TableElement(INACTIVE,7,1));
    arabidopsis.table.push_back(TableElement(INACTIVE,8,1));
    arabidopsis.table.push_back(TableElement(INACTIVE,9,1));
    arabidopsis.table.push_back(TableElement(INACTIVE,10,1));
    arabidopsis.table.push_back(TableElement(ACTIVE,1,2));
    arabidopsis.table.push_back(TableElement(INACTIVE,2,2));
    arabidopsis.table.push_back(TableElement(INACTIVE,3,2));
    arabidopsis.table.push_back(TableElement(INACTIVE,5,2));
    arabidopsis.table.push_back(TableElement(INACTIVE,8,2));
    arabidopsis.table.push_back(TableElement(INACTIVE,10,2));
    arabidopsis.table.push_back(TableElement(ACTIVE,1,3));
    arabidopsis.table.push_back(TableElement(INACTIVE,2,3));
    arabidopsis.table.push_back(TableElement(INACTIVE,5,3));
    arabidopsis.table.push_back(TableElement(INACTIVE,8,3));
    arabidopsis.table.push_back(TableElement(INACTIVE,10,3));
    arabidopsis.table.push_back(TableElement(INACTIVE,6,2));
    arabidopsis.table.push_back(TableElement(ACTIVE,1,4));
    arabidopsis.table.push_back(TableElement(ACTIVE,4,2));
    arabidopsis.table.push_back(TableElement(ACTIVE,7,2));
    arabidopsis.table.push_back(TableElement(ACTIVE,9,2));
    arabidopsis.table.push_back(TableElement(INACTIVE,2,4));
    arabidopsis.table.push_back(TableElement(INACTIVE,4,4));
    arabidopsis.table.push_back(TableElement(INACTIVE,5,4));
    arabidopsis.table.push_back(TableElement(INACTIVE,8,4));
    arabidopsis.table.push_back(TableElement(INACTIVE,10,4));
    arabidopsis.table.push_back(TableElement(INACTIVE,6,3));
    arabidopsis.table.push_back(TableElement(ACTIVE,1,5));
    arabidopsis.table.push_back(TableElement(ACTIVE,4,3));
    arabidopsis.table.push_back(TableElement(ACTIVE,7,3));
    arabidopsis.table.push_back(TableElement(ACTIVE,9,3));
    arabidopsis.table.push_back(TableElement(ACTIVE,3,3));
    arabidopsis.table.push_back(TableElement(INACTIVE,2,5));
    arabidopsis.table.push_back(TableElement(INACTIVE,4,5));
    arabidopsis.table.push_back(TableElement(INACTIVE,5,5));
    arabidopsis.table.push_back(TableElement(INACTIVE,8,5));
    arabidopsis.table.push_back(TableElement(INACTIVE,10,5));
    arabidopsis.table.push_back(TableElement(INACTIVE,6,4));
    arabidopsis.table.push_back(TableElement(ACTIVE,3,4));
    arabidopsis.table.push_back(TableElement(ACTIVE,7,4));
    arabidopsis.table.push_back(TableElement(ACTIVE,9,4));
    arabidopsis.table.push_back(TableElement(INACTIVE,6,5));
    arabidopsis.table.push_back(TableElement(ACTIVE,3,5));
    arabidopsis.table.push_back(TableElement(ACTIVE,7,5));
    arabidopsis.table.push_back(TableElement(ACTIVE,9,5));
    
    
    
    //safety check
    size_t arabidopsisTableSize = arabidopsis.geneNum * arabidopsis.timeSteps;
    size_t arabidopsisTableElementsNum = arabidopsis.table.size();
    
    if(arabidopsisTableSize != arabidopsisTableElementsNum)
    {
        std::cout << "\nERROR: missing entries in Arabidopsis timeseries table...\n\n";
        return;
    }
    
    
    
    //***********************************************************************************************************************************
    //***********************************************************************************************************************************
    //***********************************************************************************************************************************
    thcell.name = "thcell";
    thcell.timeSteps = 10;
    thcell.geneNum = 51;
    
    
    // Edges of original network
    //IRF4 		1
    //IL13 		2
    //NFAT 		3
    //MAF 		4
    //GATA3 	5
    //IL5		6
    //IFNGR		7
    //TBET		8
    //STAT1		9
    //IL7R		10
    //IFNG		11
    //STAT4		12
    //IRAK		13
    //STAT6		14
    //IL7		15
    //IL18		16
    //IL4R		17
    //IFNAR1	18
    //IFNA		19
    //IL4		20
    //SOCS1		21
    //CD80		22
    //CTLA4		23
    //SHP1		24
    //CD45		25
    //LCK		26
    //CD4		27
    //TCR		28
    //ZAP70		29
    //CD3		30
    //SLP76		31
    //VAV1		32
    //CD28		33
    //ITK		34
    //PLCPG		35
    //ANTIGEN	36
    //TNFSF4	37
    //TNFRSF4	38
    //NFKB		39
    //IKBKB		40
    //PI3K		41
    //ICOS		42
    //AKT1		43
    //COT		44
    //NIK		45
    //CD86		46
    //JAK1		47
    //JAK3		48
    //IL12R		49
    //IL12		50
    //IL18R		51
    
    thcell.edges.push_back(Edge(ACTIVATES,1,2));
    thcell.edges.push_back(Edge(ACTIVATES,3,2));
    thcell.edges.push_back(Edge(ACTIVATES,4,2));
    thcell.edges.push_back(Edge(ACTIVATES,5,2));
    thcell.edges.push_back(Edge(ACTIVATES,1,6));
    thcell.edges.push_back(Edge(ACTIVATES,3,6));
    thcell.edges.push_back(Edge(ACTIVATES,4,6));
    thcell.edges.push_back(Edge(ACTIVATES,5,6));
    thcell.edges.push_back(Edge(ACTIVATES,7,8));
    thcell.edges.push_back(Edge(INHIBITS,5,8));
    thcell.edges.push_back(Edge(ACTIVATES,9,8));
    thcell.edges.push_back(Edge(ACTIVATES,10,11));
    thcell.edges.push_back(Edge(ACTIVATES,8,11));
    thcell.edges.push_back(Edge(ACTIVATES,12,11));
    thcell.edges.push_back(Edge(ACTIVATES,9,11));
    thcell.edges.push_back(Edge(ACTIVATES,13,11));
    thcell.edges.push_back(Edge(ACTIVATES,14,4));
    thcell.edges.push_back(Edge(ACTIVATES,14,5));
    thcell.edges.push_back(Edge(INHIBITS,8,5));
    thcell.edges.push_back(Edge(ACTIVATES,15,10));
    thcell.edges.push_back(Edge(ACTIVATES,16,51));
    thcell.edges.push_back(Edge(INHIBITS,17,51));
    thcell.edges.push_back(Edge(ACTIVATES,51,13));
    thcell.edges.push_back(Edge(ACTIVATES,18,9));
    thcell.edges.push_back(Edge(ACTIVATES,7,9));
    thcell.edges.push_back(Edge(ACTIVATES,19,18));
    thcell.edges.push_back(Edge(ACTIVATES,11,7));
    thcell.edges.push_back(Edge(ACTIVATES,20,17));
    thcell.edges.push_back(Edge(INHIBITS,21,17));
    thcell.edges.push_back(Edge(ACTIVATES,1,20));
    thcell.edges.push_back(Edge(ACTIVATES,3,20));
    thcell.edges.push_back(Edge(ACTIVATES,4,20));
    thcell.edges.push_back(Edge(ACTIVATES,5,20));
    thcell.edges.push_back(Edge(ACTIVATES,22,23));
    thcell.edges.push_back(Edge(ACTIVATES,23,24));
    thcell.edges.push_back(Edge(ACTIVATES,25,26));
    thcell.edges.push_back(Edge(ACTIVATES,27,26));
    thcell.edges.push_back(Edge(ACTIVATES,28,29));
    thcell.edges.push_back(Edge(ACTIVATES,30,29));
    thcell.edges.push_back(Edge(INHIBITS,24,29));
    thcell.edges.push_back(Edge(ACTIVATES,26,29));
    thcell.edges.push_back(Edge(ACTIVATES,29,31));
    thcell.edges.push_back(Edge(ACTIVATES,26,32));
    thcell.edges.push_back(Edge(ACTIVATES,33,34));
    thcell.edges.push_back(Edge(ACTIVATES,32,34));
    thcell.edges.push_back(Edge(ACTIVATES,31,34));
    thcell.edges.push_back(Edge(ACTIVATES,34,35));
    thcell.edges.push_back(Edge(ACTIVATES,36,27));
    thcell.edges.push_back(Edge(ACTIVATES,36,28));
    thcell.edges.push_back(Edge(ACTIVATES,36,30));
    thcell.edges.push_back(Edge(ACTIVATES,36,25));
    thcell.edges.push_back(Edge(ACTIVATES,37,38));
    thcell.edges.push_back(Edge(INHIBITS,7,39));
    thcell.edges.push_back(Edge(ACTIVATES,38,39));
    thcell.edges.push_back(Edge(ACTIVATES,40,39));
    thcell.edges.push_back(Edge(ACTIVATES,14,1));
    thcell.edges.push_back(Edge(ACTIVATES,39,1));
    thcell.edges.push_back(Edge(ACTIVATES,33,3));
    thcell.edges.push_back(Edge(ACTIVATES,38,3));
    thcell.edges.push_back(Edge(ACTIVATES,35,3));
    thcell.edges.push_back(Edge(ACTIVATES,1,3));
    thcell.edges.push_back(Edge(ACTIVATES,33,41));
    thcell.edges.push_back(Edge(ACTIVATES,42,41));
    thcell.edges.push_back(Edge(ACTIVATES,41,43));
    thcell.edges.push_back(Edge(ACTIVATES,43,44));
    thcell.edges.push_back(Edge(ACTIVATES,44,45));
    thcell.edges.push_back(Edge(ACTIVATES,45,40));
    thcell.edges.push_back(Edge(ACTIVATES,46,33));
    thcell.edges.push_back(Edge(ACTIVATES,17,47));
    thcell.edges.push_back(Edge(INHIBITS,24,47));
    thcell.edges.push_back(Edge(INHIBITS,21,47));
    thcell.edges.push_back(Edge(ACTIVATES,17,48));
    thcell.edges.push_back(Edge(ACTIVATES,18,14));
    thcell.edges.push_back(Edge(ACTIVATES,47,14));
    thcell.edges.push_back(Edge(ACTIVATES,48,14));
    thcell.edges.push_back(Edge(ACTIVATES,49,12));
    thcell.edges.push_back(Edge(ACTIVATES,18,12));
    thcell.edges.push_back(Edge(INHIBITS,14,12));
    thcell.edges.push_back(Edge(ACTIVATES,7,21));
    thcell.edges.push_back(Edge(ACTIVATES,9,21));
    thcell.edges.push_back(Edge(ACTIVATES,8,21));
    thcell.edges.push_back(Edge(ACTIVATES,50,49));
    thcell.edges.push_back(Edge(ACTIVATES,37,42));
    
    
    
    
    // Added edges to corrupt the network
    thcell.addedEdges.push_back(Edge(INHIBITS,6,19));
    thcell.addedEdges.push_back(Edge(ACTIVATES,11,49));
    thcell.addedEdges.push_back(Edge(INHIBITS,20,34));
    thcell.addedEdges.push_back(Edge(ACTIVATES,41,17));
    thcell.addedEdges.push_back(Edge(ACTIVATES,39,26));
    thcell.addedEdges.push_back(Edge(ACTIVATES,23,14));
    thcell.addedEdges.push_back(Edge(ACTIVATES,43,39));
    thcell.addedEdges.push_back(Edge(ACTIVATES,34,44));
    thcell.addedEdges.push_back(Edge(INHIBITS,13,27));
    thcell.addedEdges.push_back(Edge(ACTIVATES,45,20));
    thcell.addedEdges.push_back(Edge(INHIBITS,7,43));
    thcell.addedEdges.push_back(Edge(ACTIVATES,42,12));
    thcell.addedEdges.push_back(Edge(ACTIVATES,20,38));
    thcell.addedEdges.push_back(Edge(INHIBITS,17,49));
    thcell.addedEdges.push_back(Edge(ACTIVATES,12,48));
    thcell.addedEdges.push_back(Edge(INHIBITS,48,42));
    thcell.addedEdges.push_back(Edge(INHIBITS,5,19));
    thcell.addedEdges.push_back(Edge(INHIBITS,50,36));
    thcell.addedEdges.push_back(Edge(INHIBITS,40,51));
    thcell.addedEdges.push_back(Edge(ACTIVATES,5,3));
    
    
    
    
    // Timeseries table
    thcell.table.push_back(TableElement(INACTIVE,1,1));
    thcell.table.push_back(TableElement(INACTIVE,2,1));
    thcell.table.push_back(TableElement(INACTIVE,3,1));
    thcell.table.push_back(TableElement(INACTIVE,4,1));
    thcell.table.push_back(TableElement(INACTIVE,5,1));
    thcell.table.push_back(TableElement(INACTIVE,6,1));
    thcell.table.push_back(TableElement(INACTIVE,7,1));
    thcell.table.push_back(TableElement(INACTIVE,8,1));
    thcell.table.push_back(TableElement(INACTIVE,9,1));
    thcell.table.push_back(TableElement(INACTIVE,10,1));
    thcell.table.push_back(TableElement(INACTIVE,11,1));
    thcell.table.push_back(TableElement(INACTIVE,12,1));
    thcell.table.push_back(TableElement(INACTIVE,13,1));
    thcell.table.push_back(TableElement(INACTIVE,14,1));
    thcell.table.push_back(TableElement(INACTIVE,17,1));
    thcell.table.push_back(TableElement(INACTIVE,18,1));
    thcell.table.push_back(TableElement(INACTIVE,20,1));
    thcell.table.push_back(TableElement(INACTIVE,21,1));
    thcell.table.push_back(TableElement(INACTIVE,23,1));
    thcell.table.push_back(TableElement(INACTIVE,24,1));
    thcell.table.push_back(TableElement(INACTIVE,25,1));
    thcell.table.push_back(TableElement(INACTIVE,26,1));
    thcell.table.push_back(TableElement(INACTIVE,27,1));
    thcell.table.push_back(TableElement(INACTIVE,28,1));
    thcell.table.push_back(TableElement(INACTIVE,29,1));
    thcell.table.push_back(TableElement(INACTIVE,30,1));
    thcell.table.push_back(TableElement(INACTIVE,31,1));
    thcell.table.push_back(TableElement(INACTIVE,32,1));
    thcell.table.push_back(TableElement(INACTIVE,33,1));
    thcell.table.push_back(TableElement(INACTIVE,34,1));
    thcell.table.push_back(TableElement(INACTIVE,35,1));
    thcell.table.push_back(TableElement(INACTIVE,38,1));
    thcell.table.push_back(TableElement(INACTIVE,39,1));
    thcell.table.push_back(TableElement(INACTIVE,40,1));
    thcell.table.push_back(TableElement(INACTIVE,41,1));
    thcell.table.push_back(TableElement(INACTIVE,42,1));
    thcell.table.push_back(TableElement(INACTIVE,43,1));
    thcell.table.push_back(TableElement(INACTIVE,44,1));
    thcell.table.push_back(TableElement(INACTIVE,45,1));
    thcell.table.push_back(TableElement(INACTIVE,47,1));
    thcell.table.push_back(TableElement(INACTIVE,48,1));
    thcell.table.push_back(TableElement(INACTIVE,49,1));
    thcell.table.push_back(TableElement(INACTIVE,51,1));
    thcell.table.push_back(TableElement(ACTIVE,15,1));
    thcell.table.push_back(TableElement(ACTIVE,16,1));
    thcell.table.push_back(TableElement(ACTIVE,19,1));
    thcell.table.push_back(TableElement(ACTIVE,22,1));
    thcell.table.push_back(TableElement(ACTIVE,36,1));
    thcell.table.push_back(TableElement(ACTIVE,37,1));
    thcell.table.push_back(TableElement(ACTIVE,46,1));
    thcell.table.push_back(TableElement(ACTIVE,50,1));
    thcell.table.push_back(TableElement(ACTIVE,15,2));
    thcell.table.push_back(TableElement(ACTIVE,16,2));
    thcell.table.push_back(TableElement(ACTIVE,19,2));
    thcell.table.push_back(TableElement(ACTIVE,22,2));
    thcell.table.push_back(TableElement(ACTIVE,36,2));
    thcell.table.push_back(TableElement(ACTIVE,37,2));
    thcell.table.push_back(TableElement(ACTIVE,46,2));
    thcell.table.push_back(TableElement(ACTIVE,50,2));
    thcell.table.push_back(TableElement(INACTIVE,1,2));
    thcell.table.push_back(TableElement(INACTIVE,2,2));
    thcell.table.push_back(TableElement(INACTIVE,3,2));
    thcell.table.push_back(TableElement(INACTIVE,4,2));
    thcell.table.push_back(TableElement(INACTIVE,5,2));
    thcell.table.push_back(TableElement(INACTIVE,6,2));
    thcell.table.push_back(TableElement(INACTIVE,7,2));
    thcell.table.push_back(TableElement(INACTIVE,8,2));
    thcell.table.push_back(TableElement(INACTIVE,9,2));
    thcell.table.push_back(TableElement(INACTIVE,11,2));
    thcell.table.push_back(TableElement(INACTIVE,12,2));
    thcell.table.push_back(TableElement(INACTIVE,13,2));
    thcell.table.push_back(TableElement(INACTIVE,14,2));
    thcell.table.push_back(TableElement(INACTIVE,17,2));
    thcell.table.push_back(TableElement(INACTIVE,20,2));
    thcell.table.push_back(TableElement(INACTIVE,21,2));
    thcell.table.push_back(TableElement(INACTIVE,24,2));
    thcell.table.push_back(TableElement(INACTIVE,26,2));
    thcell.table.push_back(TableElement(INACTIVE,29,2));
    thcell.table.push_back(TableElement(INACTIVE,31,2));
    thcell.table.push_back(TableElement(INACTIVE,32,2));
    thcell.table.push_back(TableElement(INACTIVE,34,2));
    thcell.table.push_back(TableElement(INACTIVE,35,2));
    thcell.table.push_back(TableElement(INACTIVE,39,2));
    thcell.table.push_back(TableElement(INACTIVE,40,2));
    thcell.table.push_back(TableElement(INACTIVE,41,2));
    thcell.table.push_back(TableElement(INACTIVE,43,2));
    thcell.table.push_back(TableElement(INACTIVE,44,2));
    thcell.table.push_back(TableElement(INACTIVE,45,2));
    thcell.table.push_back(TableElement(INACTIVE,47,2));
    thcell.table.push_back(TableElement(INACTIVE,48,2));
    thcell.table.push_back(TableElement(ACTIVE,15,3));
    thcell.table.push_back(TableElement(ACTIVE,16,3));
    thcell.table.push_back(TableElement(ACTIVE,19,3));
    thcell.table.push_back(TableElement(ACTIVE,22,3));
    thcell.table.push_back(TableElement(ACTIVE,36,3));
    thcell.table.push_back(TableElement(ACTIVE,37,3));
    thcell.table.push_back(TableElement(ACTIVE,46,3));
    thcell.table.push_back(TableElement(ACTIVE,50,3));
    thcell.table.push_back(TableElement(INACTIVE,1,3));
    thcell.table.push_back(TableElement(INACTIVE,2,3));
    thcell.table.push_back(TableElement(INACTIVE,4,3));
    thcell.table.push_back(TableElement(INACTIVE,5,3));
    thcell.table.push_back(TableElement(INACTIVE,6,3));
    thcell.table.push_back(TableElement(INACTIVE,7,3));
    thcell.table.push_back(TableElement(INACTIVE,8,3));
    thcell.table.push_back(TableElement(INACTIVE,17,3));
    thcell.table.push_back(TableElement(INACTIVE,20,3));
    thcell.table.push_back(TableElement(INACTIVE,21,3));
    thcell.table.push_back(TableElement(INACTIVE,31,3));
    thcell.table.push_back(TableElement(INACTIVE,32,3));
    thcell.table.push_back(TableElement(INACTIVE,35,3));
    thcell.table.push_back(TableElement(INACTIVE,40,3));
    thcell.table.push_back(TableElement(INACTIVE,43,3));
    thcell.table.push_back(TableElement(INACTIVE,44,3));
    thcell.table.push_back(TableElement(INACTIVE,45,3));
    thcell.table.push_back(TableElement(INACTIVE,47,3));
    thcell.table.push_back(TableElement(INACTIVE,48,3));
    thcell.table.push_back(TableElement(ACTIVE,15,4));
    thcell.table.push_back(TableElement(ACTIVE,16,4));
    thcell.table.push_back(TableElement(ACTIVE,19,4));
    thcell.table.push_back(TableElement(ACTIVE,22,4));
    thcell.table.push_back(TableElement(ACTIVE,36,4));
    thcell.table.push_back(TableElement(ACTIVE,37,4));
    thcell.table.push_back(TableElement(ACTIVE,46,4));
    thcell.table.push_back(TableElement(ACTIVE,50,4));
    thcell.table.push_back(TableElement(ACTIVE,10,2));
    thcell.table.push_back(TableElement(ACTIVE,18,2));
    thcell.table.push_back(TableElement(ACTIVE,23,2));
    thcell.table.push_back(TableElement(ACTIVE,25,2));
    thcell.table.push_back(TableElement(ACTIVE,27,2));
    thcell.table.push_back(TableElement(ACTIVE,28,2));
    thcell.table.push_back(TableElement(ACTIVE,30,2));
    thcell.table.push_back(TableElement(ACTIVE,33,2));
    thcell.table.push_back(TableElement(ACTIVE,38,2));
    thcell.table.push_back(TableElement(ACTIVE,42,2));
    thcell.table.push_back(TableElement(ACTIVE,49,2));
    thcell.table.push_back(TableElement(ACTIVE,51,2));
    thcell.table.push_back(TableElement(INACTIVE,17,4));
    thcell.table.push_back(TableElement(INACTIVE,40,4));
    thcell.table.push_back(TableElement(INACTIVE,44,4));
    thcell.table.push_back(TableElement(INACTIVE,45,4));
    thcell.table.push_back(TableElement(INACTIVE,47,4));
    thcell.table.push_back(TableElement(INACTIVE,48,4));
    thcell.table.push_back(TableElement(ACTIVE,10,3));
    thcell.table.push_back(TableElement(ACTIVE,15,5));
    thcell.table.push_back(TableElement(ACTIVE,16,5));
    thcell.table.push_back(TableElement(ACTIVE,18,3));
    thcell.table.push_back(TableElement(ACTIVE,19,5));
    thcell.table.push_back(TableElement(ACTIVE,22,5));
    thcell.table.push_back(TableElement(ACTIVE,23,3));
    thcell.table.push_back(TableElement(ACTIVE,25,3));
    thcell.table.push_back(TableElement(ACTIVE,27,3));
    thcell.table.push_back(TableElement(ACTIVE,28,3));
    thcell.table.push_back(TableElement(ACTIVE,30,3));
    thcell.table.push_back(TableElement(ACTIVE,33,3));
    thcell.table.push_back(TableElement(ACTIVE,36,5));
    thcell.table.push_back(TableElement(ACTIVE,37,5));
    thcell.table.push_back(TableElement(ACTIVE,38,3));
    thcell.table.push_back(TableElement(ACTIVE,42,3));
    thcell.table.push_back(TableElement(ACTIVE,46,5));
    thcell.table.push_back(TableElement(ACTIVE,49,3));
    thcell.table.push_back(TableElement(ACTIVE,50,5));
    thcell.table.push_back(TableElement(ACTIVE,51,3));
    thcell.table.push_back(TableElement(INACTIVE,17,5));
    thcell.table.push_back(TableElement(INACTIVE,40,5));
    thcell.table.push_back(TableElement(INACTIVE,45,5));
    thcell.table.push_back(TableElement(INACTIVE,47,5));
    thcell.table.push_back(TableElement(INACTIVE,48,5));
    thcell.table.push_back(TableElement(ACTIVE,10,4));
    thcell.table.push_back(TableElement(ACTIVE,15,6));
    thcell.table.push_back(TableElement(ACTIVE,16,6));
    thcell.table.push_back(TableElement(ACTIVE,18,4));
    thcell.table.push_back(TableElement(ACTIVE,19,6));
    thcell.table.push_back(TableElement(ACTIVE,22,6));
    thcell.table.push_back(TableElement(ACTIVE,23,4));
    thcell.table.push_back(TableElement(ACTIVE,25,4));
    thcell.table.push_back(TableElement(ACTIVE,27,4));
    thcell.table.push_back(TableElement(ACTIVE,28,4));
    thcell.table.push_back(TableElement(ACTIVE,30,4));
    thcell.table.push_back(TableElement(ACTIVE,33,4));
    thcell.table.push_back(TableElement(ACTIVE,36,6));
    thcell.table.push_back(TableElement(ACTIVE,37,6));
    thcell.table.push_back(TableElement(ACTIVE,38,4));
    thcell.table.push_back(TableElement(ACTIVE,42,4));
    thcell.table.push_back(TableElement(ACTIVE,46,6));
    thcell.table.push_back(TableElement(ACTIVE,49,4));
    thcell.table.push_back(TableElement(ACTIVE,50,6));
    thcell.table.push_back(TableElement(ACTIVE,51,4));
    thcell.table.push_back(TableElement(INACTIVE,17,6));
    thcell.table.push_back(TableElement(INACTIVE,40,6));
    thcell.table.push_back(TableElement(INACTIVE,47,6));
    thcell.table.push_back(TableElement(INACTIVE,48,6));
    thcell.table.push_back(TableElement(ACTIVE,10,5));
    thcell.table.push_back(TableElement(ACTIVE,15,7));
    thcell.table.push_back(TableElement(ACTIVE,16,7));
    thcell.table.push_back(TableElement(ACTIVE,18,5));
    thcell.table.push_back(TableElement(ACTIVE,19,7));
    thcell.table.push_back(TableElement(ACTIVE,22,7));
    thcell.table.push_back(TableElement(ACTIVE,23,5));
    thcell.table.push_back(TableElement(ACTIVE,25,5));
    thcell.table.push_back(TableElement(ACTIVE,27,5));
    thcell.table.push_back(TableElement(ACTIVE,28,5));
    thcell.table.push_back(TableElement(ACTIVE,30,5));
    thcell.table.push_back(TableElement(ACTIVE,33,5));
    thcell.table.push_back(TableElement(ACTIVE,36,7));
    thcell.table.push_back(TableElement(ACTIVE,37,7));
    thcell.table.push_back(TableElement(ACTIVE,38,5));
    thcell.table.push_back(TableElement(ACTIVE,42,5));
    thcell.table.push_back(TableElement(ACTIVE,46,7));
    thcell.table.push_back(TableElement(ACTIVE,49,5));
    thcell.table.push_back(TableElement(ACTIVE,50,7));
    thcell.table.push_back(TableElement(ACTIVE,51,5));
    thcell.table.push_back(TableElement(ACTIVE,3,3));
    thcell.table.push_back(TableElement(ACTIVE,9,3));
    thcell.table.push_back(TableElement(ACTIVE,11,3));
    thcell.table.push_back(TableElement(ACTIVE,12,3));
    thcell.table.push_back(TableElement(ACTIVE,13,3));
    thcell.table.push_back(TableElement(ACTIVE,14,3));
    thcell.table.push_back(TableElement(ACTIVE,24,3));
    thcell.table.push_back(TableElement(ACTIVE,26,3));
    thcell.table.push_back(TableElement(ACTIVE,29,3));
    thcell.table.push_back(TableElement(ACTIVE,34,3));
    thcell.table.push_back(TableElement(ACTIVE,39,3));
    thcell.table.push_back(TableElement(ACTIVE,41,3));
    thcell.table.push_back(TableElement(INACTIVE,17,7));
    thcell.table.push_back(TableElement(INACTIVE,47,7));
    thcell.table.push_back(TableElement(INACTIVE,48,7));
    thcell.table.push_back(TableElement(ACTIVE,3,4));
    thcell.table.push_back(TableElement(ACTIVE,9,4));
    thcell.table.push_back(TableElement(ACTIVE,10,6));
    thcell.table.push_back(TableElement(ACTIVE,11,4));
    thcell.table.push_back(TableElement(ACTIVE,12,4));
    thcell.table.push_back(TableElement(ACTIVE,13,4));
    thcell.table.push_back(TableElement(ACTIVE,14,4));
    thcell.table.push_back(TableElement(ACTIVE,15,8));
    thcell.table.push_back(TableElement(ACTIVE,16,8));
    thcell.table.push_back(TableElement(ACTIVE,18,6));
    thcell.table.push_back(TableElement(ACTIVE,19,8));
    thcell.table.push_back(TableElement(ACTIVE,22,8));
    thcell.table.push_back(TableElement(ACTIVE,23,6));
    thcell.table.push_back(TableElement(ACTIVE,24,4));
    thcell.table.push_back(TableElement(ACTIVE,25,6));
    thcell.table.push_back(TableElement(ACTIVE,26,4));
    thcell.table.push_back(TableElement(ACTIVE,27,6));
    thcell.table.push_back(TableElement(ACTIVE,28,6));
    thcell.table.push_back(TableElement(ACTIVE,29,4));
    thcell.table.push_back(TableElement(ACTIVE,30,6));
    thcell.table.push_back(TableElement(ACTIVE,33,6));
    thcell.table.push_back(TableElement(ACTIVE,34,4));
    thcell.table.push_back(TableElement(ACTIVE,36,8));
    thcell.table.push_back(TableElement(ACTIVE,37,8));
    thcell.table.push_back(TableElement(ACTIVE,38,6));
    thcell.table.push_back(TableElement(ACTIVE,39,4));
    thcell.table.push_back(TableElement(ACTIVE,41,4));
    thcell.table.push_back(TableElement(ACTIVE,42,6));
    thcell.table.push_back(TableElement(ACTIVE,46,8));
    thcell.table.push_back(TableElement(ACTIVE,49,6));
    thcell.table.push_back(TableElement(ACTIVE,50,8));
    thcell.table.push_back(TableElement(ACTIVE,51,6));
    thcell.table.push_back(TableElement(INACTIVE,17,8));
    thcell.table.push_back(TableElement(INACTIVE,47,8));
    thcell.table.push_back(TableElement(INACTIVE,48,8));
    thcell.table.push_back(TableElement(ACTIVE,3,5));
    thcell.table.push_back(TableElement(ACTIVE,9,5));
    thcell.table.push_back(TableElement(ACTIVE,10,7));
    thcell.table.push_back(TableElement(ACTIVE,11,5));
    thcell.table.push_back(TableElement(ACTIVE,12,5));
    thcell.table.push_back(TableElement(ACTIVE,13,5));
    thcell.table.push_back(TableElement(ACTIVE,14,5));
    thcell.table.push_back(TableElement(ACTIVE,15,9));
    thcell.table.push_back(TableElement(ACTIVE,16,9));
    thcell.table.push_back(TableElement(ACTIVE,18,7));
    thcell.table.push_back(TableElement(ACTIVE,19,9));
    thcell.table.push_back(TableElement(ACTIVE,22,9));
    thcell.table.push_back(TableElement(ACTIVE,23,7));
    thcell.table.push_back(TableElement(ACTIVE,24,5));
    thcell.table.push_back(TableElement(ACTIVE,25,7));
    thcell.table.push_back(TableElement(ACTIVE,26,5));
    thcell.table.push_back(TableElement(ACTIVE,27,7));
    thcell.table.push_back(TableElement(ACTIVE,28,7));
    thcell.table.push_back(TableElement(ACTIVE,29,5));
    thcell.table.push_back(TableElement(ACTIVE,30,7));
    thcell.table.push_back(TableElement(ACTIVE,33,7));
    thcell.table.push_back(TableElement(ACTIVE,34,5));
    thcell.table.push_back(TableElement(ACTIVE,36,9));
    thcell.table.push_back(TableElement(ACTIVE,37,9));
    thcell.table.push_back(TableElement(ACTIVE,38,7));
    thcell.table.push_back(TableElement(ACTIVE,39,5));
    thcell.table.push_back(TableElement(ACTIVE,41,5));
    thcell.table.push_back(TableElement(ACTIVE,42,7));
    thcell.table.push_back(TableElement(ACTIVE,46,9));
    thcell.table.push_back(TableElement(ACTIVE,49,7));
    thcell.table.push_back(TableElement(ACTIVE,50,9));
    thcell.table.push_back(TableElement(ACTIVE,51,7));
    thcell.table.push_back(TableElement(INACTIVE,17,9));
    thcell.table.push_back(TableElement(INACTIVE,47,9));
    thcell.table.push_back(TableElement(INACTIVE,48,9));
    thcell.table.push_back(TableElement(ACTIVE,3,6));
    thcell.table.push_back(TableElement(ACTIVE,9,6));
    thcell.table.push_back(TableElement(ACTIVE,10,8));
    thcell.table.push_back(TableElement(ACTIVE,11,6));
    thcell.table.push_back(TableElement(ACTIVE,12,6));
    thcell.table.push_back(TableElement(ACTIVE,13,6));
    thcell.table.push_back(TableElement(ACTIVE,14,6));
    thcell.table.push_back(TableElement(ACTIVE,15,10));
    thcell.table.push_back(TableElement(ACTIVE,16,10));
    thcell.table.push_back(TableElement(ACTIVE,18,8));
    thcell.table.push_back(TableElement(ACTIVE,19,10));
    thcell.table.push_back(TableElement(ACTIVE,22,10));
    thcell.table.push_back(TableElement(ACTIVE,23,8));
    thcell.table.push_back(TableElement(ACTIVE,24,6));
    thcell.table.push_back(TableElement(ACTIVE,25,8));
    thcell.table.push_back(TableElement(ACTIVE,26,6));
    thcell.table.push_back(TableElement(ACTIVE,27,8));
    thcell.table.push_back(TableElement(ACTIVE,28,8));
    thcell.table.push_back(TableElement(ACTIVE,29,6));
    thcell.table.push_back(TableElement(ACTIVE,30,8));
    thcell.table.push_back(TableElement(ACTIVE,33,8));
    thcell.table.push_back(TableElement(ACTIVE,34,6));
    thcell.table.push_back(TableElement(ACTIVE,36,10));
    thcell.table.push_back(TableElement(ACTIVE,37,10));
    thcell.table.push_back(TableElement(ACTIVE,38,8));
    thcell.table.push_back(TableElement(ACTIVE,39,6));
    thcell.table.push_back(TableElement(ACTIVE,41,6));
    thcell.table.push_back(TableElement(ACTIVE,42,8));
    thcell.table.push_back(TableElement(ACTIVE,46,10));
    thcell.table.push_back(TableElement(ACTIVE,49,8));
    thcell.table.push_back(TableElement(ACTIVE,50,10));
    thcell.table.push_back(TableElement(ACTIVE,51,8));
    thcell.table.push_back(TableElement(ACTIVE,1,4));
    thcell.table.push_back(TableElement(ACTIVE,2,4));
    thcell.table.push_back(TableElement(ACTIVE,4,4));
    thcell.table.push_back(TableElement(ACTIVE,5,4));
    thcell.table.push_back(TableElement(ACTIVE,6,4));
    thcell.table.push_back(TableElement(ACTIVE,7,4));
    thcell.table.push_back(TableElement(ACTIVE,8,4));
    thcell.table.push_back(TableElement(ACTIVE,20,4));
    thcell.table.push_back(TableElement(ACTIVE,21,4));
    thcell.table.push_back(TableElement(ACTIVE,31,4));
    thcell.table.push_back(TableElement(ACTIVE,32,4));
    thcell.table.push_back(TableElement(ACTIVE,35,4));
    thcell.table.push_back(TableElement(ACTIVE,43,4));
    thcell.table.push_back(TableElement(INACTIVE,17,10));
    thcell.table.push_back(TableElement(INACTIVE,47,10));
    thcell.table.push_back(TableElement(INACTIVE,48,10));
    thcell.table.push_back(TableElement(ACTIVE,1,5));
    thcell.table.push_back(TableElement(ACTIVE,2,5));
    thcell.table.push_back(TableElement(ACTIVE,3,7));
    thcell.table.push_back(TableElement(ACTIVE,4,5));
    thcell.table.push_back(TableElement(ACTIVE,5,5));
    thcell.table.push_back(TableElement(ACTIVE,6,5));
    thcell.table.push_back(TableElement(ACTIVE,7,5));
    thcell.table.push_back(TableElement(ACTIVE,8,5));
    thcell.table.push_back(TableElement(ACTIVE,9,7));
    thcell.table.push_back(TableElement(ACTIVE,10,9));
    thcell.table.push_back(TableElement(ACTIVE,11,7));
    thcell.table.push_back(TableElement(ACTIVE,12,7));
    thcell.table.push_back(TableElement(ACTIVE,13,7));
    thcell.table.push_back(TableElement(ACTIVE,14,7));
    thcell.table.push_back(TableElement(ACTIVE,18,9));
    thcell.table.push_back(TableElement(ACTIVE,20,5));
    thcell.table.push_back(TableElement(ACTIVE,21,5));
    thcell.table.push_back(TableElement(ACTIVE,23,9));
    thcell.table.push_back(TableElement(ACTIVE,24,7));
    thcell.table.push_back(TableElement(ACTIVE,25,9));
    thcell.table.push_back(TableElement(ACTIVE,26,7));
    thcell.table.push_back(TableElement(ACTIVE,27,9));
    thcell.table.push_back(TableElement(ACTIVE,28,9));
    thcell.table.push_back(TableElement(ACTIVE,29,7));
    thcell.table.push_back(TableElement(ACTIVE,30,9));
    thcell.table.push_back(TableElement(ACTIVE,31,5));
    thcell.table.push_back(TableElement(ACTIVE,32,5));
    thcell.table.push_back(TableElement(ACTIVE,33,9));
    thcell.table.push_back(TableElement(ACTIVE,34,7));
    thcell.table.push_back(TableElement(ACTIVE,35,5));
    thcell.table.push_back(TableElement(ACTIVE,38,9));
    thcell.table.push_back(TableElement(ACTIVE,39,7));
    thcell.table.push_back(TableElement(ACTIVE,41,7));
    thcell.table.push_back(TableElement(ACTIVE,42,9));
    thcell.table.push_back(TableElement(ACTIVE,43,5));
    thcell.table.push_back(TableElement(ACTIVE,49,9));
    thcell.table.push_back(TableElement(ACTIVE,51,9));
    thcell.table.push_back(TableElement(ACTIVE,1,6));
    thcell.table.push_back(TableElement(ACTIVE,2,6));
    thcell.table.push_back(TableElement(ACTIVE,3,8));
    thcell.table.push_back(TableElement(ACTIVE,4,6));
    thcell.table.push_back(TableElement(ACTIVE,5,6));
    thcell.table.push_back(TableElement(ACTIVE,6,6));
    thcell.table.push_back(TableElement(ACTIVE,7,6));
    thcell.table.push_back(TableElement(ACTIVE,8,6));
    thcell.table.push_back(TableElement(ACTIVE,9,8));
    thcell.table.push_back(TableElement(ACTIVE,10,10));
    thcell.table.push_back(TableElement(ACTIVE,11,8));
    thcell.table.push_back(TableElement(ACTIVE,12,8));
    thcell.table.push_back(TableElement(ACTIVE,13,8));
    thcell.table.push_back(TableElement(ACTIVE,14,8));
    thcell.table.push_back(TableElement(ACTIVE,18,10));
    thcell.table.push_back(TableElement(ACTIVE,20,6));
    thcell.table.push_back(TableElement(ACTIVE,21,6));
    thcell.table.push_back(TableElement(ACTIVE,23,10));
    thcell.table.push_back(TableElement(ACTIVE,24,8));
    thcell.table.push_back(TableElement(ACTIVE,25,10));
    thcell.table.push_back(TableElement(ACTIVE,26,8));
    thcell.table.push_back(TableElement(ACTIVE,27,10));
    thcell.table.push_back(TableElement(ACTIVE,28,10));
    thcell.table.push_back(TableElement(ACTIVE,29,8));
    thcell.table.push_back(TableElement(ACTIVE,30,10));
    thcell.table.push_back(TableElement(ACTIVE,31,6));
    thcell.table.push_back(TableElement(ACTIVE,32,6));
    thcell.table.push_back(TableElement(ACTIVE,33,10));
    thcell.table.push_back(TableElement(ACTIVE,34,8));
    thcell.table.push_back(TableElement(ACTIVE,35,6));
    thcell.table.push_back(TableElement(ACTIVE,38,10));
    thcell.table.push_back(TableElement(ACTIVE,39,8));
    thcell.table.push_back(TableElement(ACTIVE,41,8));
    thcell.table.push_back(TableElement(ACTIVE,42,10));
    thcell.table.push_back(TableElement(ACTIVE,43,6));
    thcell.table.push_back(TableElement(ACTIVE,49,10));
    thcell.table.push_back(TableElement(ACTIVE,51,10));
    thcell.table.push_back(TableElement(ACTIVE,1,7));
    thcell.table.push_back(TableElement(ACTIVE,2,7));
    thcell.table.push_back(TableElement(ACTIVE,3,9));
    thcell.table.push_back(TableElement(ACTIVE,4,7));
    thcell.table.push_back(TableElement(ACTIVE,5,7));
    thcell.table.push_back(TableElement(ACTIVE,6,7));
    thcell.table.push_back(TableElement(ACTIVE,7,7));
    thcell.table.push_back(TableElement(ACTIVE,8,7));
    thcell.table.push_back(TableElement(ACTIVE,9,9));
    thcell.table.push_back(TableElement(ACTIVE,11,9));
    thcell.table.push_back(TableElement(ACTIVE,12,9));
    thcell.table.push_back(TableElement(ACTIVE,13,9));
    thcell.table.push_back(TableElement(ACTIVE,14,9));
    thcell.table.push_back(TableElement(ACTIVE,20,7));
    thcell.table.push_back(TableElement(ACTIVE,21,7));
    thcell.table.push_back(TableElement(ACTIVE,24,9));
    thcell.table.push_back(TableElement(ACTIVE,26,9));
    thcell.table.push_back(TableElement(ACTIVE,29,9));
    thcell.table.push_back(TableElement(ACTIVE,31,7));
    thcell.table.push_back(TableElement(ACTIVE,32,7));
    thcell.table.push_back(TableElement(ACTIVE,34,9));
    thcell.table.push_back(TableElement(ACTIVE,35,7));
    thcell.table.push_back(TableElement(ACTIVE,39,9));
    thcell.table.push_back(TableElement(ACTIVE,41,9));
    thcell.table.push_back(TableElement(ACTIVE,43,7));
    thcell.table.push_back(TableElement(ACTIVE,44,5));
    thcell.table.push_back(TableElement(ACTIVE,1,8));
    thcell.table.push_back(TableElement(ACTIVE,2,8));
    thcell.table.push_back(TableElement(ACTIVE,3,10));
    thcell.table.push_back(TableElement(ACTIVE,4,8));
    thcell.table.push_back(TableElement(ACTIVE,5,8));
    thcell.table.push_back(TableElement(ACTIVE,6,8));
    thcell.table.push_back(TableElement(ACTIVE,7,8));
    thcell.table.push_back(TableElement(ACTIVE,8,8));
    thcell.table.push_back(TableElement(ACTIVE,9,10));
    thcell.table.push_back(TableElement(ACTIVE,11,10));
    thcell.table.push_back(TableElement(ACTIVE,12,10));
    thcell.table.push_back(TableElement(ACTIVE,13,10));
    thcell.table.push_back(TableElement(ACTIVE,14,10));
    thcell.table.push_back(TableElement(ACTIVE,20,8));
    thcell.table.push_back(TableElement(ACTIVE,21,8));
    thcell.table.push_back(TableElement(ACTIVE,24,10));
    thcell.table.push_back(TableElement(ACTIVE,26,10));
    thcell.table.push_back(TableElement(ACTIVE,29,10));
    thcell.table.push_back(TableElement(ACTIVE,31,8));
    thcell.table.push_back(TableElement(ACTIVE,32,8));
    thcell.table.push_back(TableElement(ACTIVE,34,10));
    thcell.table.push_back(TableElement(ACTIVE,35,8));
    thcell.table.push_back(TableElement(ACTIVE,39,10));
    thcell.table.push_back(TableElement(ACTIVE,41,10));
    thcell.table.push_back(TableElement(ACTIVE,43,8));
    thcell.table.push_back(TableElement(ACTIVE,44,6));
    thcell.table.push_back(TableElement(ACTIVE,1,9));
    thcell.table.push_back(TableElement(ACTIVE,2,9));
    thcell.table.push_back(TableElement(ACTIVE,4,9));
    thcell.table.push_back(TableElement(ACTIVE,5,9));
    thcell.table.push_back(TableElement(ACTIVE,6,9));
    thcell.table.push_back(TableElement(ACTIVE,7,9));
    thcell.table.push_back(TableElement(ACTIVE,8,9));
    thcell.table.push_back(TableElement(ACTIVE,20,9));
    thcell.table.push_back(TableElement(ACTIVE,21,9));
    thcell.table.push_back(TableElement(ACTIVE,31,9));
    thcell.table.push_back(TableElement(ACTIVE,32,9));
    thcell.table.push_back(TableElement(ACTIVE,35,9));
    thcell.table.push_back(TableElement(ACTIVE,43,9));
    thcell.table.push_back(TableElement(ACTIVE,44,7));
    thcell.table.push_back(TableElement(ACTIVE,1,10));
    thcell.table.push_back(TableElement(ACTIVE,2,10));
    thcell.table.push_back(TableElement(ACTIVE,4,10));
    thcell.table.push_back(TableElement(ACTIVE,5,10));
    thcell.table.push_back(TableElement(ACTIVE,6,10));
    thcell.table.push_back(TableElement(ACTIVE,7,10));
    thcell.table.push_back(TableElement(ACTIVE,8,10));
    thcell.table.push_back(TableElement(ACTIVE,20,10));
    thcell.table.push_back(TableElement(ACTIVE,21,10));
    thcell.table.push_back(TableElement(ACTIVE,31,10));
    thcell.table.push_back(TableElement(ACTIVE,32,10));
    thcell.table.push_back(TableElement(ACTIVE,35,10));
    thcell.table.push_back(TableElement(ACTIVE,43,10));
    thcell.table.push_back(TableElement(ACTIVE,44,8));
    thcell.table.push_back(TableElement(ACTIVE,45,6));
    thcell.table.push_back(TableElement(ACTIVE,44,9));
    thcell.table.push_back(TableElement(ACTIVE,45,7));
    thcell.table.push_back(TableElement(ACTIVE,44,10));
    thcell.table.push_back(TableElement(ACTIVE,45,8));
    thcell.table.push_back(TableElement(ACTIVE,45,9));
    thcell.table.push_back(TableElement(ACTIVE,40,7));
    thcell.table.push_back(TableElement(ACTIVE,40,8));
    thcell.table.push_back(TableElement(ACTIVE,45,10));
    thcell.table.push_back(TableElement(ACTIVE,40,9));
    thcell.table.push_back(TableElement(ACTIVE,40,10));
    
    //safety check
    size_t thcellTableSize = thcell.geneNum * thcell.timeSteps;
    size_t thcellTableElementsNum = thcell.table.size();
    
    if(thcellTableSize != thcellTableElementsNum)
    {
        std::cout << "\nERROR: missing entries in thcell timeseries table...\n\n";
        return;
    }
    
    
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    std::cout << "\n\nFINISHED LOADING NETWORKS!\n\n";
}

void LearnNetworkProperties(const std::string &repairNetwork)
{
    if(repairNetwork != "budding")
    {
        budding.LearnProperties();
        budding.PrintProperties();
    }
    else
        std::cout << "\n\n REPAIRING BUDDING NETWORK. NO LEARNING FROM BUDDING.\n\n";
    
    if(repairNetwork != "fission")
    {
        fission.LearnProperties();
        fission.PrintProperties();
    }
    else
        std::cout << "\n\n REPAIRING FISSION NETWORK. NO LEARNING FROM FISSION.\n\n";
    
    if(repairNetwork != "elegans")
    {
        elegans.LearnProperties();
        elegans.PrintProperties();
    }
    else
        std::cout << "\n\n REPAIRING ELEGANS NETWORK. NO LEARNING FROM ELEGANS.\n\n";
    
    if(repairNetwork != "mammalian")
    {
        mammalian.LearnProperties();
        mammalian.PrintProperties();
    }
    else
        std::cout << "\n\n REPAIRING MAMMALIAN NETWORK. NO LEARNING FROM MAMMALIAN.\n\n";
    
    if(repairNetwork != "arabidopsis")
    {
        arabidopsis.LearnProperties();
        arabidopsis.PrintProperties();
    }
    else
        std::cout << "\n\n REPAIRING ARABIDOPSIS NETWORK. NO LEARNING FROM ARABIDOPSIS.\n\n";
}


void AnalyzeResult(const std::string &resultFileName, const std::string &outputFileName)
{
    std::vector<Edge> originalEdges;
    
    if(resultFileName.find("budding") != std::string::npos)
        originalEdges = budding.edges;
    else if(resultFileName.find("fission") != std::string::npos)
        originalEdges = fission.edges;
    else if(resultFileName.find("elegans") != std::string::npos)
        originalEdges = elegans.edges;
    else if(resultFileName.find("mammalian") != std::string::npos)
        originalEdges = mammalian.edges;
    else if(resultFileName.find("arabidopsis") != std::string::npos)
        originalEdges = arabidopsis.edges;
    
    std::ifstream resultFile(resultFileName);
    std::ofstream outputFile(outputFileName);
    
    std::string resultLine;
    
    if(resultFile.is_open() && outputFile.is_open())
    {
        while(getline(resultFile, resultLine))
        {
            // search for Answer line
            if(resultLine.find("Answer:") != std::string::npos)
            {
                unsigned int answerNumber = 0;
                sscanf(resultLine.c_str(), "Answer: %d", &answerNumber);
                
                outputFile << "Answer: " << answerNumber << "\n";
                
                // get the actual answer set line
                getline(resultFile, resultLine);
                
                // get rid of answer sets that have "repairCost" after "activates/inhibits" (Elie's ranking case)
                std::string resultLineOnlyActivatesInhibits(resultLine, 0, resultLine.find_first_of("r")-1);
                
                resultLine = resultLineOnlyActivatesInhibits;
                
                size_t pos1 = 0;
                size_t pos2 = resultLine.find_first_of(" ");
                
                std::string restOfResultLine;
                
                std::vector<Edge> resultEdges;
                
                // stop at the edge before the last (because the end of the line is not a ' ')
                while(pos1 < pos2)
                {
                    std::string resultEdge = std::string(resultLine, pos1, pos2-pos1);
                    
                    unsigned int edgeFrom;
                    unsigned int edgeTo;
                    
                    if(resultEdge[0] == 'a')
                    {
                        sscanf(resultEdge.c_str(), "activates(%d,%d)", &edgeFrom, &edgeTo);
                        resultEdges.push_back(Edge(ACTIVATES, edgeFrom, edgeTo));
                    }
                    else
                    {
                        sscanf(resultEdge.c_str(), "inhibits(%d,%d)", &edgeFrom, &edgeTo);
                        resultEdges.push_back(Edge(INHIBITS, edgeFrom, edgeTo));
                    }
                    
                    outputFile << resultEdge << " ";
                    
                    pos1 = pos2 + 1;
                    
                    restOfResultLine = std::string(resultLine, pos1, resultLine.size());
                    pos2 += (restOfResultLine.find_first_of(" ") + 1);
                }
                
                // Don't forget last edge
                outputFile << restOfResultLine << "\n";
                
                unsigned int edgeFrom;
                unsigned int edgeTo;
                
                if(restOfResultLine[0] == 'a')
                {
                    sscanf(restOfResultLine.c_str(), "activates(%d,%d)", &edgeFrom, &edgeTo);
                    resultEdges.push_back(Edge(ACTIVATES, edgeFrom, edgeTo));
                }
                else
                {
                    sscanf(restOfResultLine.c_str(), "inhibits(%d,%d)", &edgeFrom, &edgeTo);
                    resultEdges.push_back(Edge(INHIBITS, edgeFrom, edgeTo));
                }
                
                float nbOfSimilarEdges = 0;
                
                for(size_t i = 0; i < resultEdges.size(); ++i)
                {
                    for(size_t j = 0; j < originalEdges.size(); ++j)
                    {
                        if(resultEdges[i] == originalEdges[j])
                            nbOfSimilarEdges += 1;
                    }
                }
                
                float precision = nbOfSimilarEdges / resultEdges.size();
                float recall = nbOfSimilarEdges / originalEdges.size();
                
                float f1_score = 2.0 * (precision * recall) / (precision + recall);
                
                outputFile << "\nPrecision: " << nbOfSimilarEdges << " / " << resultEdges.size() << " = " << precision << " (Nb. of edges in repaired network that are from original network)\n";
                outputFile << "Recall: " << nbOfSimilarEdges << " / " << originalEdges.size() << " = " << recall << " (Nb. of edges in original network, found in repaired network)\n";
                
                outputFile << "\nF1-score = " << f1_score << "\n";
                
                float jaccardIndex = (float)nbOfSimilarEdges / ((float)originalEdges.size() + (float)resultEdges.size() - (float)nbOfSimilarEdges);
                
                outputFile << "\nJaccard Index = " << jaccardIndex << " (Intersection of original and repaired network divided by their union)\n";
                
                outputFile << "\n";
            }
        }
        
        
        resultFile.close();
        outputFile.close();
    }
    else
    {
        std::cout << "ERROR: Unable to open result file or create output file..\n";
        return;
    }
    
    std::cout << "\n\nFINISHED ANALYZING RESULT AND CREATING FILE!\n\n";
}



std::string GetNextIntermediateName()
{
    static char character = 'A';
    static unsigned int number = 0;
    
    std::string intermediate;
    
    intermediate += character;
    intermediate += std::to_string(number);
    
    ++character;
    
    if(character > 'Z')
    {
        character = 'A';
        ++number;
    }
    
    return intermediate;
}


// ********************************************************************************************************************************************
// ********************************************************************************************************************************************
// TO DO: FINISH THIS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
// ********************************************************************************************************************************************
// ********************************************************************************************************************************************
GeneNetwork CorruptNetwork(const GeneNetwork &geneNetwork, float addedEdgesRatio, float removedEdgesRatio)
{
    GeneNetwork corruptedNetwork;
    
    if(!geneNetwork.addedEdges.empty())
    {
        std::cout << "ERROR: Trying to corrupt a network that is already corrupted!\n";
        return corruptedNetwork;
    }
    
    std::vector<Edge> edges;
    std::vector<Edge> addedEdges;
    
    for(size_t i = 0; i < geneNetwork.edges.size(); ++i)
        edges.push_back(geneNetwork.edges[i]);
    
    unsigned int nbOfEdgesToRemove = geneNetwork.edges.size() * removedEdgesRatio;
    unsigned int nbOfEdgesToAdd = geneNetwork.edges.size() * addedEdgesRatio;
    
    
    
    
    
    
    
    for(size_t i = 0; i < edges.size(); ++i)
        corruptedNetwork.edges.push_back(edges[i]);
    
    for(size_t i = 0; i < addedEdges.size(); ++i)
        corruptedNetwork.addedEdges.push_back(addedEdges[i]);
    
    return corruptedNetwork;
}



void CreateASPfile(const GeneNetwork &geneNetwork, const std::string &fileName, bool rulesOfThumb = false)
{
    std::ofstream file(fileName);
    
    if(file.is_open())
    {
        file << "% genes (nodes) of " << geneNetwork.name << " network\n";
        
        for(size_t i = 1; i <= geneNetwork.geneNum; ++i)
            file << "gene(" << i << ").\n";
        
        file << "\n% time steps\n";
        
        for(size_t i = 1; i <= geneNetwork.timeSteps; ++i)
            file << "time(" << i << ").\n";
        
        file << "\n% edges between genes\n";
        file << "% " << geneNetwork.edges.size() << " initial edges\n";
        
        size_t edgesNum = geneNetwork.edges.size();
        for(size_t currentEdge = 0; currentEdge < edgesNum; ++currentEdge)
        {
            if(geneNetwork.edges[currentEdge].type == EDGE_TYPE::ACTIVATES)
            {
                file << "edge(" << geneNetwork.edges[currentEdge].from << "," << geneNetwork.edges[currentEdge].to << ",1).\n";
            }
            
            if(geneNetwork.edges[currentEdge].type == EDGE_TYPE::INHIBITS)
            {
                file << "edge(" << geneNetwork.edges[currentEdge].from << "," << geneNetwork.edges[currentEdge].to << ",-1).\n";
            }
        }
        
        file << "\n% " << geneNetwork.addedEdges.size() << " added edges\n";
        
        if(geneNetwork.addedEdges.empty())
            std::cout << "\nWARNING: There are no added edges to corrupt the original network..\n\n";
        
        size_t addedEdgesNum = geneNetwork.addedEdges.size();
        for(size_t currentEdge = 0; currentEdge < addedEdgesNum; ++currentEdge)
        {
            if(geneNetwork.addedEdges[currentEdge].type == EDGE_TYPE::ACTIVATES)
            {
                file << "edge(" << geneNetwork.addedEdges[currentEdge].from << "," << geneNetwork.addedEdges[currentEdge].to << ",1).\n";
            }
            
            if(geneNetwork.addedEdges[currentEdge].type == EDGE_TYPE::INHIBITS)
            {
                file << "edge(" << geneNetwork.addedEdges[currentEdge].from << "," << geneNetwork.addedEdges[currentEdge].to << ",-1).\n";
            }
        }
        
        file << "\nedge(U,V) :- edge(U,V,S).\n";
        
        
        file << "\n% either add an activation edge, or an inhibition either, or nothing\n";
        file << "% this also doesn't allow the addition of both between a pair of nodes\n";
        file << "addActEdge(U,V) :- gene(U), gene(V), not edge(U,V), not addInhEdge(U,V), not nAddActEdge(U,V).\n";
        file << "nAddActEdge(U,V) :- gene(U), gene(V), not edge(U,V), not addInhEdge(U,V), not addActEdge(U,V).\n";
        
        file << "\naddInhEdge(U,V) :- gene(U), gene(V), not edge(U,V), not addActEdge(U,V), not nAddInhEdge(U,V).\n";
        file << "nAddInhEdge(U,V) :- gene(U), gene(V), not edge(U,V), not addActEdge(U,V), not addInhEdge(U,V).\n";
        
        file << "\n% either remove or don't remove existing edges\n";
        file << "removeEdge(U,V,S) :- edge(U,V,S), not nRemoveEdge(U,V,S).\n";
        file << "nRemoveEdge(U,V,S) :- edge(U,V,S), not removeEdge(U,V,S).\n";
        
        
        file << "\n% generate activates(U,V) and inhibits(U,V) to check consistency between graph and table\n";
        file << "activates(U,V) :- edge(U,V,1), not removeEdge(U,V,1).\n";
        file << "activates(U,V) :- addActEdge(U,V).\n";
        
        file << "\ninhibits(U,V) :- edge(U,V,-1), not removeEdge(U,V,-1).\n";
        file << "inhibits(U,V) :- addInhEdge(U,V).\n";
        
        file << "\n% a gene either activates or inhibits another gene\n";
        file << "% (this will make sure we don't add an activation edge and an inhibition edge at the same\n";
        file << "% time between two genes)\n";
        file << " :- activates(X,Y), inhibits(X,Y).\n";
        
        file << "\n% literals to count added edges all together\n";
        file << "% (I keep the signs so I can quickly see what is being added in answer sets)\n";
        file << "addEdge(U,V,1) :- addActEdge(U,V).\n";
        file << "addEdge(U,V,-1) :- addInhEdge(U,V).\n";
        
        file << "\n% compute cost of adding and removing edges\n";
        file << "costAdding(X) :- X = #count{addEdge(U,V,_)}.\n";
        file << "costRemoving(Y) :- Y = #count{removeEdge(U,V,_)}.\n";
        
        file << "\nrepairCost(0,Z) :- costAdding(X), costRemoving(Y), Z=X+Y.\n";
        
        
        file << "\n";
        file << "% observations at different time steps\n\n";
        file << "% timeseries table\n";
        
        size_t tableSize = geneNetwork.table.size();
        for(size_t i = 0; i < tableSize; ++i)
        {
            if(geneNetwork.table[i].type == TABLE_TYPE::ACTIVE)
            {
                file << "active(" << geneNetwork.table[i].gene << "," << geneNetwork.table[i].time << ").\n";
            }
            
            if(geneNetwork.table[i].type == TABLE_TYPE::INACTIVE)
            {
                file << "inactive(" << geneNetwork.table[i].gene << "," << geneNetwork.table[i].time << ").\n";
            }
        }
        
        
        file << "\n% activation and inhibition rules\n";
        
        file << "\n% Y receives activation at time T if X activates Y and X is active at time T\n";
        file << "receivesActivation(Y,T) :- activates(X,Y), active(X,T).\n";
        
        file << "\n% Y receives inhibition at time T if X inhibits Y and X is active at time T\n";
        file << "receivesInhibition(Y,T) :- inhibits(X,Y), active(X,T).\n";
        
        file << "\n% determine whether a gene is activated or inhibited (or none) at each time step\n";
        file << "activated(Y,T) :- receivesActivation(Y,T-1), not receivesInhibition(Y,T-1), time(T).\n";
        
        file << "\ninhibited(Y,T) :- receivesInhibition(Y,T-1), not receivesActivation(Y,T-1), time(T).\n";
        
        file << "\n% this may not be needed, but it's a sanity check\n";
        file << " :- activated(Y,T), inhibited(Y,T).\n";
        
        file << "\n% check consistency between graph and observations table\n";
        
        file << "\n% generate active and inactive for all genes based on graph propagation update\n";
        
        file << "\n% Y is inactive at time T if it was active at time T-1 and was inhibited at time T\n";
        file << "inactive(Y,T) :- active(Y,T-1), inhibited(Y,T), time(T).\n";
        
        file << "\n% Y is active at time T if it was active at time T-1 and wasn't inhibited at time T\n";
        file << "active(Y,T) :- active(Y,T-1), not inhibited(Y,T), time(T).\n";
        
        file << "\n% Y is active at time T if it was inactive at time T-1 and was activated at time T\n";
        file << "active(Y,T) :- inactive(Y,T-1), activated(Y,T), time(T).\n";
        
        file << "\n% Y is inactive at time T if it was inactive at time T-1 and wasn't activated at time T\n";
        file << "inactive(Y,T) :- inactive(Y,T-1), not activated(Y,T), time(T).\n";
        
        
        file << "\n% consistency check: make sure the data generated by the graph doesn't create conflicts\n";
        file << "% with the data in the table\n";
        file << " :- active(Y,T), inactive(Y,T).\n";
        
        if(rulesOfThumb)
        {
            
            unsigned int timeSteps = geneNetwork.timeSteps;
            unsigned int timeStepsPlus1 = timeSteps + 1;
            size_t totalEdgesNum = geneNetwork.edges.size() + geneNetwork.addedEdges.size();
            
            file << "\n\nedgeAfterRepair(U,V) :- activates(U,V).\n";
            file << "edgeAfterRepair(U,V) :- inhibits(U,V).";
            
            
            // ***********************************************************
            // RULE OF THUMB Nb. 1 => last time-step should be fixed state
            // ***********************************************************
            file << "\n\n% RULE OF THUMB Nb. 1 => last time-step should be fixed state\n";
            for(size_t i = 0; i < tableSize; ++i)
            {
                if(geneNetwork.table[i].time == geneNetwork.timeSteps)
                {
                    if(geneNetwork.table[i].type == TABLE_TYPE::ACTIVE)
                        file << "activePlus(" << geneNetwork.table[i].gene << "," << timeStepsPlus1 << ").\n";
                    
                    if(geneNetwork.table[i].type == TABLE_TYPE::INACTIVE)
                        file << "inactivePlus(" << geneNetwork.table[i].gene << "," << timeStepsPlus1 << ").\n";
                }
            }
            
            file << "\nactivated(Y," << timeStepsPlus1 << ") :- receivesActivation(Y," << timeSteps << "), not receivesInhibition(Y," << timeSteps << ").\n";
            file << "inhibited(Y," << timeStepsPlus1 << ") :- receivesInhibition(Y," << timeSteps << "), not receivesActivation(Y," << timeSteps << ").\n";
            
            file << "\nisInactivePlus(Y," << timeStepsPlus1 << ") :- active(Y," << timeSteps << "), inhibited(Y," << timeStepsPlus1 << ").\n";
            file << "isActivePlus(Y," << timeStepsPlus1 << ") :- active(Y," << timeSteps << "), not inhibited(Y," << timeStepsPlus1 << ").\n";
            file << "isActivePlus(Y," << timeStepsPlus1 << ") :- inactive(Y," << timeSteps << "), activated(Y," << timeStepsPlus1 << ").\n";
            file << "isInactivePlus(Y," << timeStepsPlus1 << ") :- inactive(Y," << timeSteps << "), not activated(Y," << timeStepsPlus1 << ").\n";
            
            file << "\n :- isActivePlus(Y," << timeStepsPlus1 << "), isInactivePlus(Y," << timeStepsPlus1 << ").\n";
            
            file << "\npenalty(1) :- activePlus(Y," << timeStepsPlus1 << "), isInactivePlus(Y," << timeStepsPlus1 << ").\n";
            file << "penalty(1) :- inactivePlus(Y," << timeStepsPlus1 << "), isActivePlus(Y," << timeStepsPlus1 << ").\n";
            
            file << "\npenalty(1) :- activePlus(Y," << timeStepsPlus1 << "), not isActivePlus(Y," << timeStepsPlus1 << ").\n";
            file << "penalty(1) :- inactivePlus(Y," << timeStepsPlus1 << "), not isInactivePlus(Y," << timeStepsPlus1 << ").\n";
            file << "penalty(1) :- not activePlus(Y," << timeStepsPlus1 << "), isActivePlus(Y," << timeStepsPlus1 << ").\n";
            file << "penalty(1) :- not inactivePlus(Y," << timeStepsPlus1 << "), isInactivePlus(Y," << timeStepsPlus1 << ").\n";
            
            file << "\nrepairCost(1," << totalEdgesNum << ") :- penalty(1).\n";
            file << "repairCost(1,0) :- not penalty(1).\n";
            
            
            // *********************************************************************************
            // RULE OF THUMB Nb. 2 => limit the number of incoming and outgoing edges for a node
            // *********************************************************************************
            file << "\n\n% RULE OF THUMB Nb. 2 => limit the number of ingoing and outgoing edges for a node\n";
            file << "kOut(C,X) :- X = #count{edgeAfterRepair(C,D)}, gene(C).\n";
            file << "kIn(C,X) :- X = #count{edgeAfterRepair(D,C)}, gene(C).\n";
            
            file << "\nkDegree(C,Z) :- kIn(C,X), kOut(C,Y), Z=X+Y.\n";
            
            file << "\nkBadGene(C) :- kDegree(C,Z), Z < 4.\n";
            file << "kBadGene(C) :- kDegree(C,Z), Z > 6.\n";
            
            file << "\nkBadGenes(X) :- X = #count{kBadGene(C)}.\n";
            
            file << "\nrepairCost(2,Y) :- kBadGenes(X), Y=X*" << (totalEdgesNum / geneNetwork.geneNum) << ".\n";
            
            
            // ***********************************************************************
            // RULE OF THUMB Nb. 3 => control the number of total edges in the network
            // ***********************************************************************
            file << "\n\n% RULE OF THUMB Nb. 3 => control the number of total edges in the network\n";
            file << "nbOfEdges(X) :- X = #count{edgeAfterRepair(C,D)}.\n";
            
            file << "\nrepairCost(3," << totalEdgesNum << "):- nbOfEdges(X), X < 21.\n";
            file << "repairCost(3," << totalEdgesNum << ") :- nbOfEdges(X), X > 25.\n";
            file << "repairCost(3,0) :- nbOfEdges(X), X >= 21, X <= 25.\n";
            
            
            // *******************************************************************************
            // RULE OF THUMB Nb. 4 => more likely interactions based on state of gene in table
            // *******************************************************************************
            file << "\n\n% RULE OF THUMB Nb. 4 => more likely interactions based on state of gene in table\n";
            
            unsigned int halfTime = geneNetwork.timeSteps / 2;
            
            file << "likelyActivator(C) :- active(C,T), T <= " << halfTime << ".\n";
            file << "likelyInhibitor(C) :- active(C,T), T > " << halfTime << ".\n";
            
            file << "\nlikelyWrongEdge(C,D) :- likelyActivator(C), inhibits(C,D), not likelyInhibitor(C), C != D.\n";
            file << "likelyWrongEdge(C,D) :- likelyInhibitor(C), activates(C,D), not likelyActivator(C), C != D.\n";
            
            file << "\nlikelyWrongEdges(X) :- X = #count{likelyWrongEdge(C,D)}.\n";
            
            file << "\nrepairCost(4,X) :- likelyWrongEdges(Y), X=Y*1.\n";
            
            
            // ***************************************************************
            // RULE OF THUMB Nb. 5 => control over the diameter of the network
            // ***************************************************************
            
            file << "\n\n% RULE OF THUMB Nb. 5 => control over the diameter of the network\n";
            file << "% first make sure every gene can be reachable\n";
            file << "link(X,Y) :- edgeAfterRepair(X,Y), X != Y.\n";
            file << "link(Y,X) :- edgeAfterRepair(X,Y), Y != X.\n";
            
            file << "\nreachable(X) :- link(1,X).\n";
            file << "reachable(Y) :- reachable(X), link(X,Y).\n";
            
            file << "\n:- gene(X), not reachable(X).\n";
            
            file << "\n% calculate diameter:\n";
            file << "% find shortest distance between every pair of vertices\n";
            file << "% diameter is the greatest value between these distances\n";
            
            for(size_t i = 0; i < 4; ++i)
            {
                file << "dist(X,Y," << (i+1) << ") :- link(X,";
                
                
                
                for(size_t j = 0; j < i; ++j)
                {
                    std::string intermediate = GetNextIntermediateName();
                    
                    file << intermediate << "), link(" << intermediate << ",";
                }
                
                file << "Y), X != Y.\n";
            }
            
            file << "\nsmallestDist(X,Y,D) :- D = #min[dist(X,Y,C)=C], dist(X,Y,Z).\n";
            
            file << "\ndiameter(D) :- D = #max[smallestDist(X,Y,C)=C].\n";
            
            file << "\nrepairCost(5," << totalEdgesNum << ") :- diameter(D), D < 3.\n";
            file << "repairCost(5," << totalEdgesNum << ") :- diameter(D), D > 4.\n";
            file << "repairCost(5,0) :- diameter(D), D >= 3, D <= 4.\n";
            
            
            // **********************************************
            // RULE OF THUMB Nb. 6 => similar dominant motifs
            // **********************************************
            
            file << "\n\n% RULE OF THUMB Nb. 6 => similar dominant motifs\n";
            file << "\nmotif3(1,X,Y,Z) :- edge(X,Y), not edge(Y,X), edge(X,Z), not edge(Z,X), not edge(Y,Z), not edge(Z,Y), X != Y, Y != Z, X != Z.\n";
            file << "motif3(2,X,Y,Z) :- not edge(X,Y), edge(Y,X), edge(X,Z), not edge(Z,X), not edge(Y,Z), not edge(Z,Y), X != Y, Y != Z, X != Z.\n";
            file << "motif3(3,X,Y,Z) :- edge(X,Y), edge(Y,X), edge(X,Z), not edge(Z,X), not edge(Y,Z), not edge(Z,Y), X != Y, Y != Z, X != Z.\n";
            file << "motif3(4,X,Y,Z) :- not edge(X,Y), not edge(Y,X), edge(X,Z), not edge(Z,X), edge(Y,Z), not edge(Z,Y), X != Y, Y != Z, X != Z.\n";
            file << "motif3(5,X,Y,Z) :- edge(X,Y), not edge(Y,X), edge(X,Z), not edge(Z,X), edge(Y,Z), not edge(Z,Y), X != Y, Y != Z, X != Z.\n";
            file << "motif3(6,X,Y,Z) :- edge(X,Y), edge(Y,X), edge(X,Z), not edge(Z,X), edge(Y,Z), not edge(Z,Y), X != Y, Y != Z, X != Z.\n";
            file << "motif3(7,X,Y,Z) :- edge(X,Y), edge(Y,X), not edge(X,Z), edge(Z,X), not edge(Y,Z), not edge(Z,Y), X != Y, Y != Z, X != Z.\n";
            file << "%motif3(8,X,Y,Z) :- edge(X,Y), edge(Y,X), edge(X,Z), edge(Z,X), not edge(Y,Z), not edge(Z,Y), X != Y, Y != Z, X != Z.\n";
            file << "%motif3(9,X,Y,Z) :- edge(X,Y), not edge(Y,X), not edge(X,Z), edge(Z,X), edge(Y,Z), not edge(Z,Y), X != Y, Y != Z, X != Z.\n";
            file << "motif3(10,X,Y,Z) :- edge(X,Y), not edge(Y,X), edge(X,Z), edge(Z,X), edge(Y,Z), not edge(Z,Y), X != Y, Y != Z, X != Z.\n";
            file << "motif3(11,X,Y,Z) :- not edge(X,Y), edge(Y,X), edge(X,Z), edge(Z,X), edge(Y,Z), not edge(Z,Y), X != Y, Y != Z, X != Z.\n";
            file << "motif3(12,X,Y,Z) :- edge(X,Y), edge(Y,X), edge(X,Z), edge(Z,X), edge(Y,Z), not edge(Z,Y), X != Y, Y != Z, X != Z.\n";
            file << "%motif3(13,X,Y,Z) :- edge(X,Y), edge(Y,X), edge(X,Z), edge(Z,X), edge(Y,Z), edge(Z,Y), X != Y, Y != Z, X != Z.\n";
            
            file << "\ndominantMotifs3(D) :- D = #count{motif3(I,X,Y,Z)}.\n";
            
            file << "\npenaltyMotifs(C) :- dominantMotifs3(Z), C=" << totalEdgesNum << "-Z.\n";
            
            file << "\nrepairCost(6,C) :- penaltyMotifs(C), C > 0.\n";
            file << "repairCost(6,0) :- penaltyMotifs(C), C <= 0.\n";
            
            
            
            
            // *********************
            // END OF RULES OF THUMB
            // *********************
            
        }
        
        
        file << "\n\n\n";
        
        if(rulesOfThumb)
        {
            file << "%totalCost(C) :- repairCost(0,C).\n";
            file << "%totalCost(C) :- repairCost(1,C).\n";
            file << "%totalCost(C) :- repairCost(2,C).\n";
            file << "%totalCost(C) :- repairCost(3,C).\n";
            file << "%totalCost(C) :- repairCost(4,C).\n";
            file << "%totalCost(C) :- repairCost(5,C).\n";
            file << "%totalCost(C) :- repairCost(6,C).\n";
        }
        
        file << "totalCost(C) :- C = #sum[repairCost(_,X)=X].\n";
        
        
        file << "\n% minimize the number of applied repairs to the network graph\n";
        file << "#minimize[totalCost(C)=C].\n";
        
        file << "\n#hide.\n";
        file << "%#show add(U,V,S).\n";
        file << "%#show remove(U,V,S).\n";
        file << "#show activates(X,Y).\n";
        file << "#show inhibits(X,Y).\n";
        file << "%#show repairCost(R,X).\n";
        file << "%#show totalCost(X).\n";
        
        
        
        file.close();
    }
    else
    {
        std::cout << "ERROR: Unable to create ASP file..\n";
        return;
    }
    
    std::cout << "\n\nFINISHED CREATING ASP FILE!\n\n";
}



void StatisticalApproachWithSignificance(const std::string &randomRepairsFileName, const std::string &outputFileName)
{
    std::vector<Edge> originalEdges;
    
    if(randomRepairsFileName.find("budding") != std::string::npos)
        originalEdges = budding.edges;
    else if(randomRepairsFileName.find("fission") != std::string::npos)
        originalEdges = fission.edges;
    else if(randomRepairsFileName.find("elegans") != std::string::npos)
        originalEdges = elegans.edges;
    else if(randomRepairsFileName.find("mammalian") != std::string::npos)
        originalEdges = mammalian.edges;
    else if(randomRepairsFileName.find("arabidopsis") != std::string::npos)
        originalEdges = arabidopsis.edges;
    
    std::ifstream randomRepairsFile(randomRepairsFileName);
    std::ofstream outputFile(outputFileName);
    
    std::string repairLine;
    
    std::vector<Repair> repairs;
    
    unsigned int repairNumber = 0;
    
    if(randomRepairsFile.is_open() && outputFile.is_open())
    {
        while(getline(randomRepairsFile, repairLine))
        {
            // search for Answer line
            if(repairLine.find("Answer:") != std::string::npos)
            {
                outputFile << "\nRepair: " << repairNumber++ << "\n\n";
                
                // get the actual answer set line
                getline(randomRepairsFile, repairLine);
                
                long startHere = repairLine.find_first_of("r");
                
                std::string repairEdgesString(repairLine, 0, startHere-1);
                std::string repairCostString(repairLine, startHere, repairLine.size()-1);
                
                unsigned int rule0 = 0;
                unsigned int rule1 = 0;
                unsigned int rule2 = 0;
                unsigned int rule3 = 0;
                unsigned int rule4 = 0;
                unsigned int rule5 = 0;
                unsigned int rule6 = 0;
                
                // save rule costs for statistical approach
                sscanf(repairCostString.c_str(), "repairCost(0,%d) repairCost(1,%d) repairCost(2,%d) repairCost(3,%d) repairCost(4,%d) repairCost(5,%d) repairCost(6,%d)", &rule0, &rule1, &rule2, &rule3, &rule4, &rule5, &rule6);
                
                repairs.push_back(Repair(rule0, rule1, rule2, rule3, rule4, rule5, rule6));
                
                // evaluate each repair, then pick the best one based on statistical approach
                // I decided to evaluate all repairs in case we need to do some comparisons...
                size_t pos1 = 0;
                size_t pos2 = repairEdgesString.find_first_of(" ");
                
                std::string restOfRepairEdgesString;
                
                std::vector<Edge> resultEdges;
                
                // stop at the edge before the last (because the end of the line is not a ' ')
                while(pos1 < pos2)
                {
                    std::string resultEdge = std::string(repairEdgesString, pos1, pos2-pos1);
                    
                    unsigned int edgeFrom;
                    unsigned int edgeTo;
                    
                    if(resultEdge[0] == 'a')
                    {
                        sscanf(resultEdge.c_str(), "activates(%d,%d)", &edgeFrom, &edgeTo);
                        resultEdges.push_back(Edge(ACTIVATES, edgeFrom, edgeTo));
                    }
                    else
                    {
                        sscanf(resultEdge.c_str(), "inhibits(%d,%d)", &edgeFrom, &edgeTo);
                        resultEdges.push_back(Edge(INHIBITS, edgeFrom, edgeTo));
                    }
                    
                    outputFile << resultEdge << " ";
                    
                    pos1 = pos2 + 1;
                    
                    restOfRepairEdgesString = std::string(repairEdgesString, pos1, repairEdgesString.size());
                    pos2 += (restOfRepairEdgesString.find_first_of(" ") + 1);
                }
                
                // Don't forget last edge
                outputFile << restOfRepairEdgesString << "\n";
                
                unsigned int edgeFrom;
                unsigned int edgeTo;
                
                if(restOfRepairEdgesString[0] == 'a')
                {
                    sscanf(restOfRepairEdgesString.c_str(), "activates(%d,%d)", &edgeFrom, &edgeTo);
                    resultEdges.push_back(Edge(ACTIVATES, edgeFrom, edgeTo));
                }
                else
                {
                    sscanf(restOfRepairEdgesString.c_str(), "inhibits(%d,%d)", &edgeFrom, &edgeTo);
                    resultEdges.push_back(Edge(INHIBITS, edgeFrom, edgeTo));
                }
                
                float nbOfSimilarEdges = 0;
                
                for(size_t i = 0; i < resultEdges.size(); ++i)
                {
                    for(size_t j = 0; j < originalEdges.size(); ++j)
                    {
                        if(resultEdges[i] == originalEdges[j])
                            nbOfSimilarEdges += 1;
                    }
                }
                
                float precision = nbOfSimilarEdges / resultEdges.size();
                float recall = nbOfSimilarEdges / originalEdges.size();
                
                float f1_score = 2.0 * (precision * recall) / (precision + recall);
                
                outputFile << "\nPrecision: " << nbOfSimilarEdges << " / " << resultEdges.size() << " = " << precision << " (Nb. of edges in repaired network that are from original network)\n";
                outputFile << "Recall: " << nbOfSimilarEdges << " / " << originalEdges.size() << " = " << recall << " (Nb. of edges in original network, found in repaired network)\n";
                
                outputFile << "\nF1-score = " << f1_score << "\n";
                
                float jaccardIndex = (float)nbOfSimilarEdges / ((float)originalEdges.size() + (float)resultEdges.size() - (float)nbOfSimilarEdges);
                
                outputFile << "\nJaccard Index = " << jaccardIndex << " (Intersection of original and repaired network divided by their union)\n";
                
                outputFile << "\n";
            }
        }
        
        // Pick the best repair based on statistical approach
        
        float averages[7] = {0};
        
        size_t repairsNum = repairs.size();
        for(size_t i = 0; i < repairsNum; ++i)
        {
            averages[0] += repairs[i].ruleViolations[0];
            averages[1] += repairs[i].ruleViolations[1];
            averages[2] += repairs[i].ruleViolations[2];
            averages[3] += repairs[i].ruleViolations[3];
            averages[4] += repairs[i].ruleViolations[4];
            averages[5] += repairs[i].ruleViolations[5];
            averages[6] += repairs[i].ruleViolations[6];
        }
        
        for(size_t i = 0; i < 7; ++i)
            averages[i] /= (float)repairsNum;
        
        
        float variances[7] = {0};
        
        for(size_t i = 0; i < repairsNum; ++i)
        {
            float value0 = (float)(repairs[i].ruleViolations[0]) - averages[0];
            variances[0] += (value0 * value0);
            
            float value1 = (float)(repairs[i].ruleViolations[1]) - averages[1];
            variances[1] += (value1 * value1);
            
            float value2 = (float)(repairs[i].ruleViolations[2]) - averages[2];
            variances[2] += (value2 * value2);
            
            float value3 = (float)(repairs[i].ruleViolations[3]) - averages[3];
            variances[3] += (value3 * value3);
            
            float value4 = (float)(repairs[i].ruleViolations[4]) - averages[4];
            variances[4] += (value4 * value4);
            
            float value5 = (float)(repairs[i].ruleViolations[5]) - averages[5];
            variances[5] += (value5 * value5);
            
            float value6 = (float)(repairs[i].ruleViolations[6]) - averages[6];
            variances[6] += (value6 * value6);
        }
        
        for(size_t i = 0; i < 7; ++i)
            variances[i] /= ((float)repairsNum - 1.0);
        
        float standardDeviations[7] = {0};
        
        for(size_t i = 0; i < 7; ++i)
            standardDeviations[i] = sqrtf(variances[i]);
        
        
        for(size_t i = 0; i < repairsNum; ++i)
        {
            for(size_t j = 0; j < 7; ++j)
            {
                if(standardDeviations[j] <= 0.000001f)
                {
                    repairs[i].zScores[j] = 0;
                    continue;
                }
                
                repairs[i].zScores[j] = ((float)(repairs[i].ruleViolations[j]) - averages[j]) / standardDeviations[j];
                repairs[i].totalZScore += repairs[i].zScores[j];
            }
        }
        
        // pick repair with smallest z-score as best repair
        float smallestZScore = 999999.0f;
        size_t bestRepair = 0;
        
        for(size_t i = 0; i < repairsNum; ++i)
        {
            if(repairs[i].totalZScore < smallestZScore)
            {
                smallestZScore = repairs[i].totalZScore;
                bestRepair = i;
            }
        }
        
        outputFile << "\n\n\nBEST REPAIR:\n";
        outputFile << "============\n\n";
        
        outputFile << "Repair: " << bestRepair;
        
        outputFile << "\n\nCost of rules: ";
        
        for(int i = 0; i < 7; ++i)
            outputFile << repairs[bestRepair].ruleViolations[i] << " ";
        
        outputFile << "\n\n\n";
        
        randomRepairsFile.close();
        outputFile.close();
    }
    else
    {
        std::cout << "ERROR: Unable to open random repairs file or create output file..\n";
        return;
    }
    
    std::cout << "\n\nFINISHED STATISTICAL APPROACH AND CREATED OUTPUT FILE!\n\n";
}

void StatisticalApproachWithSignificance2(const std::string &randomRepairsFileName, const std::string &outputFileName)
{
    std::vector<Edge> originalEdges;
    
    if(randomRepairsFileName.find("budding") != std::string::npos)
        originalEdges = budding.edges;
    else if(randomRepairsFileName.find("fission") != std::string::npos)
        originalEdges = fission.edges;
    else if(randomRepairsFileName.find("elegans") != std::string::npos)
        originalEdges = elegans.edges;
    else if(randomRepairsFileName.find("mammalian") != std::string::npos)
        originalEdges = mammalian.edges;
    else if(randomRepairsFileName.find("arabidopsis") != std::string::npos)
        originalEdges = arabidopsis.edges;
    
    std::ifstream randomRepairsFile(randomRepairsFileName);
    std::ofstream outputFile(outputFileName);
    
    std::string repairLine;
    
    std::vector<Repair> repairs;
    
    if(randomRepairsFile.is_open() && outputFile.is_open())
    {
        while(getline(randomRepairsFile, repairLine))
        {
            // search for Answer line
            if(repairLine.find("Answer:") != std::string::npos)
            {
                // get the actual answer set line
                getline(randomRepairsFile, repairLine);
                
                long startHere = repairLine.find_first_of("r");
                
                std::string repairEdgesString(repairLine, 0, startHere-1);
                std::string repairCostString(repairLine, startHere, repairLine.size()-1);
                
                unsigned int rule0 = 0;
                unsigned int rule1 = 0;
                unsigned int rule2 = 0;
                unsigned int rule3 = 0;
                unsigned int rule4 = 0;
                unsigned int rule5 = 0;
                unsigned int rule6 = 0;
                
                // save rule costs for statistical approach
                sscanf(repairCostString.c_str(), "repairCost(0,%d) repairCost(1,%d) repairCost(2,%d) repairCost(3,%d) repairCost(4,%d) repairCost(5,%d) repairCost(6,%d)", &rule0, &rule1, &rule2, &rule3, &rule4, &rule5, &rule6);
                
                repairs.push_back(Repair(rule0, rule1, rule2, rule3, rule4, rule5, rule6));
                
            }
        }
        
        // Calculate averages and standard deviations
        
        float averages[7] = {0};
        
        size_t repairsNum = repairs.size();
        for(size_t i = 0; i < repairsNum; ++i)
        {
            averages[0] += repairs[i].ruleViolations[0];
            averages[1] += repairs[i].ruleViolations[1];
            averages[2] += repairs[i].ruleViolations[2];
            averages[3] += repairs[i].ruleViolations[3];
            averages[4] += repairs[i].ruleViolations[4];
            averages[5] += repairs[i].ruleViolations[5];
            averages[6] += repairs[i].ruleViolations[6];
        }
        
        for(size_t i = 0; i < 7; ++i)
            averages[i] /= (float)repairsNum;
        
        
        float variances[7] = {0};
        
        for(size_t i = 0; i < repairsNum; ++i)
        {
            float value0 = (float)(repairs[i].ruleViolations[0]) - averages[0];
            variances[0] += (value0 * value0);
            
            float value1 = (float)(repairs[i].ruleViolations[1]) - averages[1];
            variances[1] += (value1 * value1);
            
            float value2 = (float)(repairs[i].ruleViolations[2]) - averages[2];
            variances[2] += (value2 * value2);
            
            float value3 = (float)(repairs[i].ruleViolations[3]) - averages[3];
            variances[3] += (value3 * value3);
            
            float value4 = (float)(repairs[i].ruleViolations[4]) - averages[4];
            variances[4] += (value4 * value4);
            
            float value5 = (float)(repairs[i].ruleViolations[5]) - averages[5];
            variances[5] += (value5 * value5);
            
            float value6 = (float)(repairs[i].ruleViolations[6]) - averages[6];
            variances[6] += (value6 * value6);
        }
        
        for(size_t i = 0; i < 7; ++i)
            variances[i] /= ((float)repairsNum - 1.0);
        
        float standardDeviations[7] = {0};
        
        for(size_t i = 0; i < 7; ++i)
            standardDeviations[i] = sqrtf(variances[i]);
        
        
        outputFile << "\n\nAVERAGES:\n";
        outputFile << "=========\n\n";
        
        outputFile << "rule 0: " << averages[0] << "\n";
        outputFile << "rule 1: " << averages[1] << "\n";
        outputFile << "rule 2: " << averages[2] << "\n";
        outputFile << "rule 3: " << averages[3] << "\n";
        outputFile << "rule 4: " << averages[4] << "\n";
        outputFile << "rule 5: " << averages[5] << "\n";
        outputFile << "rule 6: " << averages[6] << "\n";
        
        outputFile << "\n\nSTANDARD DEVIATIONS:\n";
        outputFile << "====================\n\n";
        
        outputFile << "rule 0: " << standardDeviations[0] << "\n";
        outputFile << "rule 1: " << standardDeviations[1] << "\n";
        outputFile << "rule 2: " << standardDeviations[2] << "\n";
        outputFile << "rule 3: " << standardDeviations[3] << "\n";
        outputFile << "rule 4: " << standardDeviations[4] << "\n";
        outputFile << "rule 5: " << standardDeviations[5] << "\n";
        outputFile << "rule 6: " << standardDeviations[6] << "\n";
        
        outputFile << "\n\nASP CODE:\n";
        outputFile << "=========\n\n";
        
        std::vector<std::string> zScoreNames;
        
        for(int i = 0; i < 7; ++i)
        {
            if((averages[i] < 0.001) || (standardDeviations[i] < 0.1))
                continue;
            
            outputFile << "cost" << i << "(X) :- repairCost(" << i << ",C), X=1000*C.\n";
            outputFile << "diff" << i << "(D) :- cost" << i << "(X), D=X-" << (int)(averages[i]*1000) << ".\n";
            outputFile << "zScore" << i << "(Z) :- diff" << i << "(D), Z=D/" << (int)(standardDeviations[i]*10) << ".\n\n";
            
            std::string zScoreName;
            
            zScoreName += "zScore";
            zScoreName += std::to_string(i);
            
            zScoreNames.push_back(zScoreName);
        }
        
        if(zScoreNames.size() < 2)
        {
            std::cout << "ERROR: There is only 1 z-score greater than 0..\n";
            return;
        }
        
        
        outputFile << "totalCost0(X) :- " << zScoreNames[0] << "(A), " << zScoreNames[1] << "(B), X=A+B.\n";
        
        size_t zScoreNum = zScoreNames.size();
        size_t c;
        
        for(c = 2; c < zScoreNum; ++c)
        {
            outputFile << "totalCost" << (c - 1) << "(X) :- totalCost" << (c - 2) << "(A), " << zScoreNames[c] << "(B), X=A+B.\n";
        }
        
        outputFile << "\ntotalCost(X) :- totalCost" << (c - 2) << "(X).\n\n";
        
        outputFile << "#minimize[totalCost(C)=C].\n\n";
        
        outputFile << "#hide.\n";
        outputFile << "#show activates(X,Y).\n";
        outputFile << "#show inhibits(X,Y).\n";
        
        
        randomRepairsFile.close();
        outputFile.close();
    }
    else
    {
        std::cout << "ERROR: Unable to open random repairs file or create output file..\n";
        return;
    }
    
    std::cout << "\n\nFINISHED STATISTICAL APPROACH AND CREATED OUTPUT FILE!\n\n";
}


void ElieRanking(const std::string &aspFileName, const std::string &outputFileName)
{
    std::string bestRepairFileName;
    
    if(aspFileName.find("budding") != std::string::npos)
        bestRepairFileName = "bestRepair_buddingElieRanking.txt";
    else if(aspFileName.find("fission") != std::string::npos)
        bestRepairFileName = "bestRepair_fissionElieRanking.txt";
    else if(aspFileName.find("elegans") != std::string::npos)
        bestRepairFileName = "bestRepair_elegansElieRanking.txt";
    else if(aspFileName.find("mammalian") != std::string::npos)
        bestRepairFileName = "bestRepair_mammalianElieRanking.txt";
    else if(aspFileName.find("arabidopsis") != std::string::npos)
        bestRepairFileName = "bestRepair_arabidopsisElieRanking.txt";
    
    bool timeLimitReached = false;
    unsigned int changeCounter = 0;
    
    while(!timeLimitReached)
    {
        int values[7] = {0};
        
        // run ASP solvers to get the next "better" answer set repair
        std::string commandString;
        
        commandString += "./gringo ";
        commandString += aspFileName;
        commandString += " | clasp --time-limit=10 > repair.txt";
        
        std::system(commandString.c_str());
        
        std::string repairLine;
        
        // check if time limit was reached before the solvers could find a better repair
        std::ifstream repairFileCopyCheck("repair.txt");
        
        if(repairFileCopyCheck.is_open())
        {
            while(getline(repairFileCopyCheck, repairLine))
            {
                // if time limit was reached, stop everything and exit
                if(repairLine.find("UNKNOWN") != std::string::npos)
                {
                    std::cout << "\n\nTime limit reached. bestRepair.txt file contains the best repair found. Exiting..\n\n\n";
                    repairFileCopyCheck.close();
                    timeLimitReached = true;
                    
                    std::cout << "\n\nFINISHED ELIE'S RANKING APPROACH AND CREATED OUTPUT FILE!\n\n";
                    
                    // analyze the result we get
                    AnalyzeResult(bestRepairFileName, outputFileName);
                    
                    return;
                }
            }
        }
        else
        {
            std::cout << "ERROR: Unable to open repair file..\n";
            return;
        }
        
        // if time limit was not reached, a repair was found, make a copy of it (this is the best repair so far)
        std::ifstream repairFileCopy("repair.txt");
        std::ofstream bestRepairFile(bestRepairFileName);
        
        if(repairFileCopy.is_open() && bestRepairFile.is_open())
        {
            while(getline(repairFileCopy, repairLine))
            {
                bestRepairFile << repairLine << std::endl;
                std::cout << repairLine << std::endl;
            }
            
            repairFileCopy.close();
            repairFileCopyCheck.close();
            bestRepairFile.close();
        }
        else
        {
            std::cout << "ERROR: Unable to open repair file..\n";
            return;
        }
        
        std::ifstream repairFile("repair.txt");
        
        // read repair to get the penalty values of each rule
        if(repairFile.is_open())
        {
            while(getline(repairFile, repairLine))
            {
                // search for Answer line
                if(repairLine.find("Answer:") != std::string::npos)
                {
                    // get the actual answer set line
                    getline(repairFile, repairLine);
                    
                    long startHere = repairLine.find_first_of("r");
                    
                    std::string repairEdgesString(repairLine, 0, startHere-1);
                    std::string repairCostString(repairLine, startHere, repairLine.size()-1);
                    
                    unsigned int rule0 = 0;
                    unsigned int rule1 = 0;
                    unsigned int rule2 = 0;
                    unsigned int rule3 = 0;
                    unsigned int rule4 = 0;
                    unsigned int rule5 = 0;
                    unsigned int rule6 = 0;
                    
                    // save rule penalties
                    sscanf(repairCostString.c_str(), "repairCost(0,%d) repairCost(1,%d) repairCost(2,%d) repairCost(3,%d) repairCost(4,%d) repairCost(5,%d) repairCost(6,%d)", &rule0, &rule1, &rule2, &rule3, &rule4, &rule5, &rule6);
                    
                    values[0] = rule0;
                    values[1] = rule1;
                    values[2] = rule2;
                    values[3] = rule3;
                    values[4] = rule4;
                    values[5] = rule5;
                    values[6] = rule6;
                }
            }
            
            repairFile.close();
        }
        else
        {
            std::cout << "ERROR: Unable to open repair file..\n";
            return;
        }
        
        // write new ASP file containing updated rule penalties
        std::ifstream aspFile(aspFileName);
        std::ofstream aspFileNew("newASP.txt");
        
        std::string aspLine;
        
        if(aspFile.is_open() && aspFileNew.is_open())
        {
            while(getline(aspFile, aspLine))
            {
                if(aspLine.find("#hide") == std::string::npos)
                    aspFileNew << aspLine << std::endl;
                else
                    break;
            }
            
            aspFileNew << "change" << changeCounter << "(0,-1) :- repairCost(0,X), X < " << values[0] << ".\n";
            aspFileNew << "change" << changeCounter << "(1,-1) :- repairCost(1,X), X < " << values[1] << ".\n";
            aspFileNew << "change" << changeCounter << "(2,-1) :- repairCost(2,X), X < " << values[2] << ".\n";
            aspFileNew << "change" << changeCounter << "(3,-1) :- repairCost(3,X), X < " << values[3] << ".\n";
            aspFileNew << "change" << changeCounter << "(4,-1) :- repairCost(4,X), X < " << values[4] << ".\n";
            aspFileNew << "change" << changeCounter << "(5,-1) :- repairCost(5,X), X < " << values[5] << ".\n";
            aspFileNew << "change" << changeCounter << "(6,-1) :- repairCost(6,X), X < " << values[6] << ".\n\n";
            
            aspFileNew << "change" << changeCounter << "(0,1) :- repairCost(0,X), X > " << values[0] << ".\n";
            aspFileNew << "change" << changeCounter << "(1,1) :- repairCost(1,X), X > " << values[1] << ".\n";
            aspFileNew << "change" << changeCounter << "(2,1) :- repairCost(2,X), X > " << values[2] << ".\n";
            aspFileNew << "change" << changeCounter << "(3,1) :- repairCost(3,X), X > " << values[3] << ".\n";
            aspFileNew << "change" << changeCounter << "(4,1) :- repairCost(4,X), X > " << values[4] << ".\n";
            aspFileNew << "change" << changeCounter << "(5,1) :- repairCost(5,X), X > " << values[5] << ".\n";
            aspFileNew << "change" << changeCounter << "(6,1) :- repairCost(6,X), X > " << values[6] << ".\n\n";
            
            aspFileNew << "totalChange" << changeCounter << "(C) :- C = #sum[change" << changeCounter << "(X,Y)=Y].\n\n";
            
            aspFileNew << " :- totalChange" << changeCounter << "(C), C >= 0.\n\n";
            
            aspFileNew << "#hide.\n";
            aspFileNew << "#show repairCost(X,Y).\n";
            aspFileNew << "#show activates(X,Y).\n";
            aspFileNew << "#show inhibits(X,Y).\n";
            
            ++changeCounter;
            
            
            // make the new ASP file as input file for next iteration
            std::rename(aspFileName.c_str(), "oldASP.txt");
            std::rename("newASP.txt", aspFileName.c_str());
            
            aspFile.close();
            aspFileNew.close();
        }
        else
        {
            std::cout << "ERROR: Unable to open ASP file..\n";
            return;
        }
    }
    
    // analyze the result we get
    AnalyzeResult(bestRepairFileName, outputFileName);
    
    std::cout << "\n\nSOLVERS DID NOT TIME OUT... NO REPAIR WAS FOUND WITHIN THE TIME LIMIT!\n\n";
    std::cout << "\n\nFINISHED ELIE'S RANKING APPROACH AND CREATED OUTPUT FILE!\n\n";
}

void FindAverages(const std::string &fileName)
{
    std::ifstream file(fileName);
    
    std::string line;
    
    float numbers[10] = {0};
    float averages[10] = {0};
    
    
    if(file.is_open())
    {
        while(getline(file, line))
        {
            sscanf(line.c_str(), "coordinates {(1, %f)(2, %f)(3, %f)(4, %f)(5, %f)(6, %f)(7, %f)(8, %f)(9, %f)(10, %f)};",
                   &(numbers[0]), &(numbers[1]), &(numbers[2]), &(numbers[3]), &(numbers[4]), &(numbers[5]),
                   &(numbers[6]), &(numbers[7]), &(numbers[8]), &(numbers[9]));
            
            for(int i = 0; i < 10; ++i)
                averages[i] += numbers[i];
        }
        
        for(int i = 0; i < 10; ++i)
            averages[i] /= 5.0f;
        
        for(int i = 0; i < 10; ++i)
            std::cout << "  " << numbers[i];
        
        std::cout << std::endl;
        
        for(int i = 0; i < 10; ++i)
            std::cout << std::setprecision(2) << "(" << i+1 << ", " << averages[i] << ")";
        
        std::cout << std::endl;
        
    }
}



int main(int argc, const char * argv[])
{
    //    LoadNetworks(NOT_CORRUPTED);
    //    LearnNetworkProperties("budding");
    
    
    
    
    
    
    //    LoadNetworks(CORRUPTED);
    //
    //    CreateASPfile(budding, "budding.txt", true);
    //    CreateASPfile(budding, "buddingNoRules.txt", false);
    //
    //    CreateASPfile(fission, "fission.txt", true);
    //    CreateASPfile(fission, "fissionNoRules.txt", false);
    //
    //    CreateASPfile(elegans, "elegans.txt", true);
    //    CreateASPfile(elegans, "elegansNoRules.txt", false);
    //
    //    CreateASPfile(mammalian, "mammalian.txt", true);
    //    CreateASPfile(mammalian, "mammalianNoRules.txt", false);
    //
    //    CreateASPfile(arabidopsis, "arabidopsis.txt", true);
    //    CreateASPfile(arabidopsis, "arabidopsisNoRules.txt", false);
    
    
    
    
    
    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! WARNING: answers set should contain ONLY "activates(X,Y)" and "inhibits(X,Y)" predicates
    
    LoadNetworks(NOT_CORRUPTED);
    
    
    //ElieRanking("mammalianElieRanking.txt", "FINALRESULT_mammalianElieRanking.txt");
    
    //AnalyzeResult("bestRepair_mammalianElieRanking.txt","FINALRESULT_mammalianVioImp.txt");
    
    //AnalyzeResult("result_mammalianNoRules_cons.txt","FINALRESULT_mammalianNoRules.txt");
    //AnalyzeResult("result_mammalian_cons.txt","FINALRESULT_mammalian.txt");
    //AnalyzeResult("result_mammalianLeximin_cons.txt","FINALRESULT_mammalianLeximin.txt");
    //AnalyzeResult("result_mammalianLeximax_cons.txt","FINALRESULT_mammalianLeximax.txt");
    //AnalyzeResult("result_mammalianStatistical_cons.txt","FINALRESULT_mammalianStatistical.txt");
    //AnalyzeResult("result_mammalianStatisticalMinimal_cons.txt","FINALRESULT_mammalianStatisticalMinimal.txt");
    //AnalyzeResult("result_mammalianStatisticalThreeMinimal_cons.txt","FINALRESULT_mammalianStatisticalThreeMinimal.txt");
    
    
    //
    //    AnalyzeResult("result_buddingNoRules_80_20.txt","FINALRESULT_buddingNoRules_80_20.txt");
    //    AnalyzeResult("result_budding_80_20.txt","FINALRESULT_budding_80_20.txt");
    //    AnalyzeResult("result_buddingLeximin_80_20.txt","FINALRESULT_buddingLeximin_80_20.txt");
    //    AnalyzeResult("result_buddingLeximax_80_20.txt","FINALRESULT_buddingLeximax_80_20.txt");
    //    AnalyzeResult("result_buddingStatistical_80_20.txt","FINALRESULT_buddingStatistical_80_20.txt");
    //    AnalyzeResult("result_buddingStatisticalMinimal_80_20.txt","FINALRESULT_buddingStatisticalMinimal_80_20.txt");
    //    AnalyzeResult("result_buddingStatisticalThreeMinimal_80_20.txt","FINALRESULT_buddingStatisticalThreeMinimal_80_20.txt");
    
    //    AnalyzeResult("result_fissionNoRules_80_20.txt","FINALRESULT_fissionNoRules_80_20.txt");
    //    AnalyzeResult("result_fission_80_20.txt","FINALRESULT_fission_80_20.txt");
    //    AnalyzeResult("result_fissionLeximin_80_20.txt","FINALRESULT_fissionLeximin_80_20.txt");
    //    AnalyzeResult("result_fissionLeximax_80_20.txt","FINALRESULT_fissionLeximax_80_20.txt");
    //    AnalyzeResult("result_fissionStatistical_80_20.txt","FINALRESULT_fissionStatistical_80_20.txt");
    //    AnalyzeResult("result_fissionStatisticalMinimal_80_20.txt","FINALRESULT_fissionStatisticalMinimal_80_20.txt");
    //    AnalyzeResult("result_fissionStatisticalThreeMinimal_80_20.txt","FINALRESULT_fissionStatisticalThreeMinimal_80_20.txt");
    
    //    AnalyzeResult("result_elegansNoRules_80_20.txt","FINALRESULT_elegansNoRules_80_20.txt");
    //    AnalyzeResult("result_elegans_80_20.txt","FINALRESULT_elegans_80_20.txt");
    //    AnalyzeResult("result_elegansLeximin_80_20.txt","FINALRESULT_elegansLeximin_80_20.txt");
    //    AnalyzeResult("result_elegansLeximax_80_20.txt","FINALRESULT_elegansLeximax_80_20.txt");
    //    AnalyzeResult("result_elegansStatistical_80_20.txt","FINALRESULT_elegansStatistical_80_20.txt");
    //    AnalyzeResult("result_elegansStatisticalMinimal_80_20.txt","FINALRESULT_elegansStatisticalMinimal_80_20.txt");
    //    AnalyzeResult("result_elegansStatisticalThreeMinimal_80_20.txt","FINALRESULT_elegansStatisticalThreeMinimal_80_20.txt");
    
    //    AnalyzeResult("result_mammalianNoRules_80_20.txt","FINALRESULT_mammalianNoRules_80_20.txt");
    //    AnalyzeResult("result_mammalian_80_20.txt","FINALRESULT_mammalian_80_20.txt");
    //    AnalyzeResult("result_mammalianLeximin_80_20.txt","FINALRESULT_mammalianLeximin_80_20.txt");
    //    AnalyzeResult("result_mammalianLeximax_80_20.txt","FINALRESULT_mammalianLeximax_80_20.txt");
    //    AnalyzeResult("result_mammalianStatistical_80_20.txt","FINALRESULT_mammalianStatistical_80_20.txt");
    //    AnalyzeResult("result_mammalianStatisticalMinimal_80_20.txt","FINALRESULT_mammalianStatisticalMinimal_80_20.txt");
    //    AnalyzeResult("result_mammalianStatisticalThreeMinimal_80_20.txt","FINALRESULT_mammalianStatisticalThreeMinimal_80_20.txt");
    
    
    //    AnalyzeResult("result_arabidopsisNoRules_80_20.txt","FINALRESULT_arabidopsisNoRules_80_20.txt");
    //    AnalyzeResult("result_arabidopsis_80_20.txt","FINALRESULT_arabidopsis_80_20.txt");
    //    AnalyzeResult("result_arabidopsisLeximin_80_20.txt","FINALRESULT_arabidopsisLeximin_80_20.txt");
    //    AnalyzeResult("result_arabidopsisLeximax_80_20.txt","FINALRESULT_arabidopsisLeximax_80_20.txt");
    //    AnalyzeResult("result_arabidopsisStatistical_80_20.txt","FINALRESULT_arabidopsisStatistical_80_20.txt");
    //    AnalyzeResult("result_arabidopsisStatisticalMinimal_80_20.txt","FINALRESULT_arabidopsisStatisticalMinimal_80_20.txt");
    //    AnalyzeResult("result_arabidopsisStatisticalThreeMinimal_80_20.txt","FINALRESULT_arabidopsisStatisticalThreeMinimal_80_20.txt");
    
    
    
    
    
    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! WARNING: answer sets should contain ONLY "activates(X,Y)", "inhibits(X,Y)" and "repairCost(X,Y)" predicates
    
    //    LoadNetworks(NOT_CORRUPTED);
    
    //    StatisticalApproachWithSignificance("dispersedRepairs_budding_80_20.txt", "FINALRESULT_buddingStatistical_80_20.txt");
    //    StatisticalApproachWithSignificance("dispersedRepairs_fission_80_20.txt", "FINALRESULT_fissionStatistical_80_20.txt");
    //    StatisticalApproachWithSignificance("dispersedRepairs_elegans_80_20.txt", "FINALRESULT_elegansStatistical_80_20.txt");
    //    StatisticalApproachWithSignificance("dispersedRepairs_mammalian_80_20.txt", "FINALRESULT_mammalianStatistical_80_20.txt");
    //    StatisticalApproachWithSignificance("dispersedRepairs_arabidopsis_80_20.txt", "FINALRESULT_arabidopsisStatistical_80_20.txt");
    
    
    
    
    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! WARNING: answer sets should contain ONLY "repairCost(X,Y)" predicates
    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! WARNING: I use the ASP code in the output of this to create ASP file that will find best repairs
    
    //    LoadNetworks(NOT_CORRUPTED);
    
    //    StatisticalApproachWithSignificance2("randomRepairs_budding_80_20.txt", "ASP_buddingStatistical.txt");
    //    StatisticalApproachWithSignificance2("randomRepairs_fission_80_20.txt", "ASP_fissionStatistical.txt");
    //    StatisticalApproachWithSignificance2("randomRepairs_elegans_80_20.txt", "ASP_elegansStatistical.txt");
    //    StatisticalApproachWithSignificance2("randomRepairs_mammalian_80_20.txt", "ASP_mammalianStatistical.txt");
    //    StatisticalApproachWithSignificance2("randomRepairs_arabidopsis_80_20.txt", "ASP_arabidopsisStatistical.txt");
    
    
    //    StatisticalApproachWithSignificance2("minimalRepairs_budding_80_20.txt", "ASP_buddingStatistical_80_20.txt");
    //    StatisticalApproachWithSignificance2("minimalRepairs_fission_80_20.txt", "ASP_fissionStatistical_80_20.txt");
    //    StatisticalApproachWithSignificance2("minimalRepairs_elegans_80_20.txt", "ASP_elegansStatistical_80_20.txt");
    //    StatisticalApproachWithSignificance2("minimalRepairs_mammalian_80_20.txt", "ASP_mammalianStatistical_80_20.txt");
    //    StatisticalApproachWithSignificance2("minimalRepairs_arabidopsis_80_20.txt", "ASP_arabidopsisStatistical_80_20.txt");
    
    
    //    StatisticalApproachWithSignificance2("threeMinimalRepairs_budding_80_20.txt", "ASP_buddingStatistical_80_20.txt");
    //    StatisticalApproachWithSignificance2("threeMinimalRepairs_fission_80_20.txt", "ASP_fissionStatistical_80_20.txt");
    //    StatisticalApproachWithSignificance2("threeMinimalRepairs_elegans_80_20.txt", "ASP_elegansStatistical_80_20.txt");
    //    StatisticalApproachWithSignificance2("threeMinimalRepairs_mammalian_80_20.txt", "ASP_mammalianStatistical_80_20.txt");
    //    StatisticalApproachWithSignificance2("threeMinimalRepairs_arabidopsis_80_20.txt", "ASP_arabidopsisStatistical_80_20.txt");
    
    
    
    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! WARNING: Straight forward. This automatically calls gringo and clasp and give the final result from here.
    // !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! WARNING: Should have gringo in the debug folder and clasp on the machine path (/usr/bin).
    
    //LoadNetworks(NOT_CORRUPTED);
    
    
    //ElieRanking("buddingElieRanking.txt", "FINALRESULT_buddingElieRanking_80_20.txt");
    //ElieRanking("fissionElieRanking.txt", "FINALRESULT_fissionElieRanking_80_20.txt");
    //ElieRanking("elegansElieRanking.txt", "FINALRESULT_elegansElieRanking_80_20.txt");
    //ElieRanking("mammalianElieRanking.txt", "FINALRESULT_mammalianElieRanking_80_20.txt");
    //ElieRanking("arabidopsisElieRanking.txt", "FINALRESULT_arabidopsisElieRanking_80_20.txt");
    
    //ElieRanking("buddingElieRanking.txt", "FINALRESULT_buddingElieRanking_cons.txt");
    //ElieRanking("fissionElieRanking.txt", "FINALRESULT_fissionElieRanking_cons.txt");
    //ElieRanking("elegansElieRanking.txt", "FINALRESULT_elegansElieRanking_cons.txt");
    //ElieRanking("mammalianElieRanking.txt", "FINALRESULT_mammalianElieRanking_cons.txt");
    //ElieRanking("arabidopsisElieRanking.txt", "FINALRESULT_arabidopsisElieRanking_cons.txt");
    
    
    
    
    
    return 0;
}







