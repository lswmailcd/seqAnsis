#include"Common.h"
#include"AppException.h"
#include"GlobalSpace.h"
#include"PhylogeneticTree.h"

namespace SeqAnsis
{

	CPhylogeneticTree::CPhylogeneticTree() :m_root(NULL), m_nptr(NULL), m_ptrs(NULL), m_lptr(NULL), m_olptr(NULL)
{
}

CPhylogeneticTree::~CPhylogeneticTree()
{
	clearTree(m_root);
}

bool CPhylogeneticTree::readTree(const std::string& strTreeFileName, const std::vector<CSequence>& Seqs)
{
	if (m_root)
	{
		clearTree(m_root);
	}

	m_file.open(strTreeFileName);
	if (!m_file)
	{
		return false;
	}

	nnodes = 0;
	ntotal = 0;
	nSeq = 0;

	m_nSeq = (int)Seqs.size();

	if (m_nptr)	 SAFE_DELETE_ARRAY(m_nptr);
	m_nptr = new TreeNode*[3 * m_nSeq];
	if (m_ptrs)	 SAFE_DELETE_ARRAY(m_ptrs);
	m_ptrs = new TreeNode*[3 * m_nSeq];
	if (m_lptr)	 SAFE_DELETE_ARRAY(m_lptr);
	m_lptr = new TreeNode*[m_nSeq];
	if (m_olptr)	 SAFE_DELETE_ARRAY(m_olptr);
	m_olptr = new TreeNode*[m_nSeq];

	m_root = avail();
	setInfo(m_root, NULL, 0, std::string(""), 0.0);

	createTree(m_root, NULL, &m_file);

	m_file.close();

	if (nSeq != m_nSeq)
	{
		char ch[100];
		sprintf_s(ch, "tree not compatible with alignment. %d sequences in alignment and %d in tree", m_nSeq, nSeq);
		CGlobalSpace::m_sEventLog.writeEvent(ch);
		return false;
	}
	// If the tree is unrooted, reroot the tree - ie. minimise the difference
	// between the mean root->leaf distances for the left and right branches of
	// the tree.     
#if 0
	if (clustalw::userParameters->getDistanceTree() == false)
	{
		if (rootedTree == false)
		{
			clustalw::utilityObject->error("input tree is unrooted and has no distances.\nCannot align sequences");
			return 0;
		}
	}

	if (rootedTree == false)
	{
		root = reRoot(seqTree, lastSeq - firstSeq + 1);
	}
	else
	{
		root = seqTree;
	}
#endif

	// calculate the 'order' of each node.
	int nameLength;
	std::string nameSeq;
	orderNodes();

	if (m_nSeq >= 2)
	{
		// If there are more than three sequences....
		// assign the sequence nodes (in the same order as in the alignment file)
		for (int i = 0; i < m_nSeq; i++)
		{
			nameLength = Seqs[i].getName().length();
			nameSeq = Seqs[i].getName();
			std::string name1 = "";
			std::string name2 = "";
			if (nameLength > MAXNAMES)
			{
				char ch[200];
				sprintf_s(ch, "name %s is too long for PHYLIP tree format", nameSeq.c_str());
				CGlobalSpace::m_sEventLog.writeEvent(ch);
			}

			for (int k = 0; k < nameLength && k < MAXNAMES; k++)
			{
				char c = nameSeq[k];
				if ((c > 0x40) && (c < 0x5b))
				{
					c = c | 0x20;
				}
				if (c == ' ')
				{
					c = '_';
				}
				name2 += c;
			}
			bool found = false;

			for (int j = 0; j < m_nSeq; j++)
			{
				name1 = "";
				for (int k = 0; k < (int)m_lptr[j]->name.length() && k < MAXNAMES; k++)
				{
					char c = m_lptr[j]->name[k];
					if ((c > 0x40) && (c < 0x5b))
					{
						c = c | 0x20;
					}

					name1 += c;
				}

				if (name1.compare(name2) == 0)
				{
					m_olptr[i] = m_lptr[j];
					found = true;
				}
			}

			if (found == false)
			{
				char ch[100];
				sprintf_s(ch, "tree not compatible with alignment!  %s not found\n", name2.c_str());
				CGlobalSpace::m_sEventLog.writeEvent(ch);
				return false;
			}
		}

	}
	return true;
}

TreeNode* CPhylogeneticTree::insertRoot(TreeNode* p, float diff)
{
	TreeNode *newp, *prev, *q, *t;
	float dist, prevDist, td;

	newp = avail();

	if (p->parent == NULL)
	{
		// AW bug 94: question remains if access here should be handled differently
		char ch[100];
		sprintf_s(ch, "INTERNAL ERROR: Tree::insertRoot: TreeNode p->parent is NULL");
		CGlobalSpace::m_sEventLog.writeEvent(ch);
		throw CAppException(DEF_EXCEPTION_NULL_POINTER, DEF_EXCEPTION_LEVEL_EXIT_TRY_BLOCK, __EXCEPTION_SITE__, ch);
	}
	
	t = p->parent;
	prevDist = t->dist;

	p->parent = newp;

	dist = p->dist;

	p->dist = diff / 2;

	if (p->dist < 0.0)
	{
		p->dist = 0.0;
	}

	if (p->dist > dist)
	{
		p->dist = dist;
	}

	t->dist = dist - p->dist;

	newp->left = t;
	newp->right = p;
	newp->parent = NULL;
	newp->dist = 0.0;
	newp->leaf = NODE;
	if (t->left == p)
	{
		t->left = t->parent;
	}
	else
	{
		t->right = t->parent;
	}
	prev = t;
	q = t->parent;

	t->parent = newp;

	while (q != NULL)
	{
		if (q->left == prev)
		{
			q->left = q->parent;
			q->parent = prev;
			td = q->dist;
			q->dist = prevDist;
			prevDist = td;
			prev = q;
			q = q->left;
		}
		else
		{
			q->right = q->parent;
			q->parent = prev;
			td = q->dist;
			q->dist = prevDist;
			prevDist = td;
			prev = q;
			q = q->right;
		}
	}

		/*
		* remove the old root node
		*/
	q = prev;
	if (q->left == NULL)
	{
		dist = q->dist;
		q = q->right;
		q->dist += dist;
		q->parent = prev->parent;
		if (prev->parent->left == prev)
		{
			prev->parent->left = q;
		}
		else
		{
			prev->parent->right = q;
		}

		prev->right = NULL;
	}
	else
	{
		dist = q->dist;
		q = q->left;
		q->dist += dist;
		q->parent = prev->parent;

		if (prev->parent->left == prev)
		{
			prev->parent->left = q;
		}
		else
		{
			prev->parent->right = q;
		}

		prev->left = NULL;
	}

	return newp;
}

void CPhylogeneticTree::calcSeqWeights(int nSeqNum, std::vector<float>& SeqWeight)
{
	if ((int)SeqWeight.size() < nSeqNum)
    {
		SeqWeight.resize(nSeqNum);
    }
    
    int i;
    int temp;
	float *weight, sum;
    /*
     * If there are more than three sequences....
     */

	if (nSeqNum >= 2)
    {
        /*
         * Calculate sequence weights based on Phylip tree.
         */

		weight = new float[nSeqNum];
        
		for (i = 0; i < nSeqNum; ++i)
        {
            weight[i] = calcWeight(i);
        }

        /*
         * Normalise the weights(0,1), such that the sum of the weights = 1.0f
         */

        sum = 0;
		for (i = 0; i < nSeqNum; ++i)
        {
            sum += weight[i];
        }

        if (sum == 0)
        {
			for (i = 0; i < nSeqNum; ++i)
            {
                weight[i] = 1.0f;
            }
			sum = (float)nSeqNum;
        }

		for (i = 0; i < nSeqNum; ++i)
        {
			SeqWeight[i] = weight[i] / sum;
        }        
		SAFE_DELETE_ARRAY(weight);
    }

    else
    {
		for (i = 0; i < nSeqNum; ++i)
        {
			SeqWeight[i] = 1.0f / nSeqNum;
        }
    }
}

/** *************************************************************************
 * Private functions!!!!!!!!!!!!!!!                                         *
 ****************************************************************************/
 
/**
 * 
 * @param ptree 
 * @param parent 
 * @param file 
 */
void CPhylogeneticTree::createTree(TreeNode* ptree, TreeNode* parent, std::ifstream* file)
{
    TreeNode* p;

    int i, type;
    float dist;
    std::string name;

    
    // is this a node or a leaf ?
    skipSpace(file);
	m_charFromFile = file->get();
	if (m_charFromFile == '(')
    {
        // this must be a node....
        type = NODE;
        name = "";
        m_ptrs[ntotal] = m_nptr[nnodes] = ptree;
        nnodes++;
        ntotal++;

        createNode(ptree, parent);

        p = ptree->left;
        createTree(p, ptree, file);

		if (m_charFromFile == ',')
        {
            p = ptree->right;
            createTree(p, ptree, file);
			if (m_charFromFile == ',')
            {
                ptree = insertNode(ptree);
				m_ptrs[ntotal] = m_nptr[nnodes] = ptree;
                nnodes++;
                ntotal++;
                p = ptree->right;
                createTree(p, ptree, file);
                rootedTree = false;
            }
        }

        skipSpace(file);
		m_charFromFile = file->get();
    }
    
    // ...otherwise, this is a leaf
    
    else
    {
        type = LEAF;
		m_ptrs[ntotal++] = m_lptr[nSeq++] = ptree;
        // get the sequence name
        name = "";
		name += m_charFromFile;
		m_charFromFile = file->get();
        
        i = 1;
		while ((m_charFromFile != ':') && (m_charFromFile != ',') && (m_charFromFile != ')'))
        {
            if (i < MAXNAMES)
            {
				name += m_charFromFile;
                i++;
            }
			m_charFromFile = file->get();
        }

		if (m_charFromFile != ':')
        {
            dist = 0.0;
        }
    }
    
    // get the distance information
    
    dist = 0.0;
	if (m_charFromFile == ':')
    {
        skipSpace(file);
        (*file) >> dist;
        skipSpace(file);
		m_charFromFile = file->get();
    }
    setInfo(ptree, parent, type, name, dist);
}

/**
 * 
 * @param pptr 
 * @param parent 
 */
void CPhylogeneticTree::createNode(TreeNode* pptr, TreeNode* parent)
{
    TreeNode* t;
    pptr->parent = parent;
    t = avail();
    pptr->left = t;
    t = avail();
    pptr->right = t;
}

/**
 * 
 * @param pptr 
 * @return 
 */
TreeNode* CPhylogeneticTree::insertNode(TreeNode* pptr)
{
	TreeNode* newnode = avail();
    createNode(newnode, pptr->parent);

    newnode->left = pptr;
    pptr->parent = newnode;

    setInfo(newnode, pptr->parent, NODE, "", 0.0);

    return newnode;
}

/**
 * 
 * @param p 
 */
void CPhylogeneticTree::clearTree(TreeNode* p)
{
	if (!p)  return;

    clearTreeNodes(p);
    delete [] m_nptr;
	m_nptr = NULL;
	delete[] m_ptrs;
	m_ptrs = NULL;
	delete[] m_lptr;
	m_lptr = NULL;
	delete[] m_olptr;
	m_olptr = NULL;
}

/**
 * 
 * @param p 
 */
void CPhylogeneticTree::clearTreeNodes(TreeNode* p)
{
    if (p == NULL)
    {
        p = m_root;
    }
    if (p->left != NULL)
    {
        clearTreeNodes(p->left);
    }
    if (p->right != NULL)
    {
        clearTreeNodes(p->right);
    }
    p->left = NULL;
    p->right = NULL;
    
    delete p;
    p  = NULL;
}

/**
 * 
 * @param ptree 
 * @param nseqs 
 * @return 
 */
TreeNode* CPhylogeneticTree::reRoot(TreeNode* ptree, int nseqs)
{
    TreeNode *p, *rootNode, *rootPtr;
    float diff, minDiff = 0.0, minDepth = 1.0, maxDist;
    int i;
    bool first = true;

    // find the difference between the means of leaf->node
    // distances on the left and on the right of each node
    rootPtr = ptree;
    for (i = 0; i < ntotal; i++)
    {
		p = m_ptrs[i];
        if (p->parent == NULL)
        {
            /* AW Bug 94: p->parent must be chosen as rootNode
               (non-optimized executables (-O0) never do), otherwise
               insertRoot fails.
               Is p->parent == NULL valid at all?
               Simplest thing for now is to continue here. Tree code
               needs serious dismantling anyway. See debugPrintAllNodes
            */

            continue;
            //diff = calcRootMean(p, &maxDist);
        }
        else
        {
            diff = calcMean(p, &maxDist, nseqs);
        }

        if ((diff == 0) || ((diff > 0) && (diff < 2 *p->dist)))
        {
            if ((maxDist < minDepth) || (first == true))
            {
                first = false;
                rootPtr = p;
                minDepth = maxDist;
                minDiff = diff;
            }
        }

    }

    
    // insert a new node as the ancestor of the node which produces the shallowest
    // tree.
    /* AW Bug 94: could also be prevented here */
    if (rootPtr == ptree)
    {
        minDiff = rootPtr->left->dist + rootPtr->right->dist;
        rootPtr = rootPtr->right;
    }
    rootNode = insertRoot(rootPtr, minDiff);

    diff = calcRootMean(rootNode, &maxDist);

    return rootNode;
}

/**
 * 
 * @param root 
 * @param maxDist 
 * @return 
 */
float CPhylogeneticTree::calcRootMean(TreeNode* root, float *maxDist)
{
    float dist, leftSum = 0.0, rightSum = 0.0, leftMean, rightMean, diff;
    TreeNode* p;
    int i;
    int numLeft, numRight;
    int direction;
    
    // for each leaf, determine whether the leaf is left or right of the root.

    dist = (*maxDist) = 0;
    numLeft = numRight = 0;
	for (i = 0; i < m_nSeq; i++)
    {
        p = m_lptr[i];
        dist = 0.0;
        while (p->parent != root)
        {
            dist += p->dist;
            p = p->parent;
        }

        if (p == root->left)
        {
            direction = LEFT;
        }
        else
        {
            direction = RIGHT;
        }

        dist += p->dist;

        if (direction == LEFT)
        {
            leftSum += dist;
            numLeft++;
        }
        else
        {
            rightSum += dist;
            numRight++;
        }
        if (dist > (*maxDist))
        {
            *maxDist = dist;
        }
    }

    leftMean = leftSum / numLeft;
    rightMean = rightSum / numRight;

    diff = leftMean - rightMean;
    return diff;
}

/**
 * 
 * @param nptr 
 * @param maxDist 
 * @param nSeqs 
 * @return 
 */
float CPhylogeneticTree::calcMean(TreeNode* nptr, float *maxDist, int nSeqs)
{
    float dist, leftSum = 0.0, rightSum = 0.0, leftMean, rightMean, diff;
    TreeNode* p;  
    TreeNode** pathToRoot;
    float *distToNode;
    int depth = 0, i, j, n = 0;
    int numLeft, numRight;
    int direction;
    bool found;

    pathToRoot = new TreeNode*[nSeqs];
    distToNode = new float[nSeqs];
    // determine all nodes between the selected node and the root;
  
    depth = 0;
    (*maxDist) = dist = 0.0;
    numLeft = numRight = 0;
    p = nptr;
    while (p != NULL)
    {
        pathToRoot[depth] = p;
        dist += p->dist;
        distToNode[depth] = dist;
        p = p->parent;
        depth++;
    }

    // for each leaf, determine whether the leaf is left or right of the node.
    // (RIGHT = descendant, LEFT = not descendant)
 
	for (i = 0; i < m_nSeq; i++)
    {
        p = m_lptr[i];
        if (p == nptr)
        {
            direction = RIGHT;
            dist = 0.0;
        }
        else
        {
            direction = LEFT;
            dist = 0.0;
            
            // find the common ancestor. 
            
            found = false;
            n = 0;
            while ((found == false) && (p->parent != NULL))
            {
                for (j = 0; j < depth; j++)
                if (p->parent == pathToRoot[j])
                {
                    found = true;
                    n = j;
                }

                dist += p->dist;
                p = p->parent;
            }
            if (p == nptr)
            {
                direction = RIGHT;
            }
        }

        if (direction == LEFT)
        {
            leftSum += dist;
            leftSum += distToNode[n - 1];
            numLeft++;
        }
        else
        {
            rightSum += dist;
            numRight++;
        }

        if (dist > (*maxDist))
        {
            *maxDist = dist;
        }
    }

    delete [] distToNode;
    distToNode = NULL;
    delete [] pathToRoot;
    pathToRoot = NULL;
    
    leftMean = leftSum / numLeft;
    rightMean = rightSum / numRight;

    diff = leftMean - rightMean;
    return (diff);
}

/**
 * 
 */
void CPhylogeneticTree::orderNodes()
{
    int i;
    TreeNode* p;

	for (i = 0; i < m_nSeq; i++)
    {
        p = m_lptr[i];
        while (p != NULL)
        {
            p->order++;
            p = p->parent;
        }
    }
}

/**
 * 
 * @param leaf 
 * @return 
 */
float CPhylogeneticTree::calcWeight(int leaf)
{
    float weight = 0.0;
	TreeNode*p = m_olptr[leaf];
    
	while (p->parent != NULL)
    {
        weight += p->dist / p->order;
        p = p->parent;
    }

    //weight *= 100.0;

    return weight;
}

/**
 * skipSpace is used to skip all the spaces at the begining of a file. The next read
 * will give a character other than a space.
 * @param file 
 */
void CPhylogeneticTree::skipSpace(std::ifstream* file)
{
    char c;

    do
    {
        c = file->get();
    }
    while (isspace(c));

    file->putback(c);
}

/**
 * 
 * @param p 
 * @param groups 
 * @param n 
 */
void CPhylogeneticTree::markGroup1(TreeNode* p, int *groups, int n)
{
    int i;

    for (i = 0; i < n; i++)
    {
        if (m_olptr[i] == p)
        {
            groups[i] = 1;
        }
        else
        {
            groups[i] = 0;
        }
    }
}

/**
 * 
 * @param p 
 * @param groups 
 * @param n 
 */
void CPhylogeneticTree::markGroup2(TreeNode* p, int *groups, int n)
{
    int i;

    for (i = 0; i < n; i++)
    {
        if (m_olptr[i] == p)
        {
            groups[i] = 2;
        }
        else if (groups[i] != 0)
        {
            groups[i] = 1;
        }
    }
}

/**
 * 
 * @return 
 */
TreeNode* CPhylogeneticTree::avail()
{
    TreeNode* p;
    p = new TreeNode; 
    p->left = NULL;
    p->right = NULL;
    p->parent = NULL;
    p->dist = 0.0;
    p->leaf = 0;
    p->order = 0;
    return (p);
}

void CPhylogeneticTree::setInfo(TreeNode* p, TreeNode* parent, int pleaf, const std::string& pname, float pdist)
{
    p->parent = parent;
    p->leaf = pleaf;
    p->dist = pdist;
    p->order = 0;
    p->name = pname;
    if (p->leaf == true)
    {
        p->left = NULL;
        p->right = NULL;
    }
}
          
}
