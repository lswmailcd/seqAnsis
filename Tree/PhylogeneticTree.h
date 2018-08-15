#pragma once

#include <vector>
#include <string>
#include <iostream>
#include <fstream>
#include <memory>
#include"Sequence.h"

namespace SeqAnsis
{
	const int MAXNAMES = 100;
	const int INT_SCALE_FACTOR = 100;
	const int NODE = 0;
	const int LEAF = 1;
	const int LEFT = 1;
	const int RIGHT = 2;

struct TreeNode
{
// phylogenetic tree structure
	struct TreeNode *left;
	struct TreeNode *right;
	struct TreeNode *parent;
	float dist;
	int leaf;
	int order;
	std::string name;
};

class CPhylogeneticTree
{
    public:
		CPhylogeneticTree();
		~CPhylogeneticTree();
        void calcSeqWeights(int nSeqNum, std::vector<float>& SeqWeight);
        void clearTree(TreeNode* p);
		bool readTree(const std::string& strTreeFileName, const std::vector<CSequence>& Seqs);
    private:
        void createTree(TreeNode* ptree, TreeNode* parent, std::ifstream* file);
        void createNode(TreeNode* pptr, TreeNode* parent);
        TreeNode* insertNode(TreeNode* pptr);
        void clearTreeNodes(TreeNode* p);
        TreeNode* reRoot(TreeNode* ptree, int nseqs);
        TreeNode* insertRoot(TreeNode* p, float diff);
        float calcRootMean(TreeNode* root, float *maxDist);
        float calcMean(TreeNode* nptr, float *maxDist, int nSeqs);
        void orderNodes();
        float calcWeight(int leaf);
        void skipSpace(std::ifstream* file);
        void markGroup1(TreeNode* p, int *groups, int n);
        void markGroup2(TreeNode* p, int *groups, int n);
        TreeNode* avail();
        void setInfo(TreeNode* p, TreeNode* parent, int pleaf, const std::string& pname, float pdist);
            
        char m_charFromFile;
        std::ifstream m_file;
        TreeNode** m_lptr;
        TreeNode** m_olptr;
        TreeNode** m_nptr;
        TreeNode** m_ptrs;
        int nnodes;
        int ntotal;
		int nSeq;
        bool rootedTree;
        TreeNode* seqTree;
        TreeNode* m_root;
		int m_nSeq;
        int* groups;
        int numSets;
};

}