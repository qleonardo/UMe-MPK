#ifndef DC_H
#define DC_H

#include <stdint.h>
#include <string>
#include <vector>
#include <cstdio>

// D&C tree structure
struct DCNode 
{
	int *level_offset;
	int nbParts, firstRow, lastRow;
    struct DCNode **son, *iso;
	bool isLeaf;

	~DCNode()
	{
		if (level_offset != NULL)
			delete[] level_offset;
		if (son != NULL)
			delete[] son;
		if (iso != NULL)
			delete iso;
	}
};

// D&C arguments structure
struct DCArgs 
{
    int firstRow, lastRow;
};

typedef struct 
{
    int *index, *value;
} index_t;

using namespace std;

class DC
{
	private:
		int *RowPerm, *RowRev;
		int PARTSIZE, LEVEL, Recursion;
	
		int DC_partitioning (int **Row2Row, int *nRowPerRow, int *RowValue, int **local_Row2Row, int *local_nRowPerRow, int firstRow, int lastRow, int *RowPart, bool flag);

		void DC_create (DCNode *tree, int **Row2Row, int *RowValue, int *nRowPerRow, int globalNbRow, int firstRow, int lastRow, int recursion);

		void DC_create_level(DCNode *tree, int **Row2Row, int *nRowPerRow, int *RowValue, int firstRow, int lastRow);
									  
		// Apply local element permutation to global element permutation
		void merge_permutations (int *perm, int *localPerm, int globalNbRow, int localNbRow, int firstRow, int lastRow);

	public:
	
		DCNode *treeRoot;

		void DC_creation (int *c2c, int *index, int globalNbRow);

		void DC_traversal (void(*userFctPtrA) (char **, DCArgs *), void(*userFctPtrB) (char **, 
                                    DCArgs *), char ***userArgs, int total);

		void DC_traversal_level (void(*userFctPtrA) (char **, DCArgs *), 
									void(*userFctPtrB) (char **, DCArgs *), 
									char ***userArgs, int total);

		int DC_get_blocks(int m, int **row_offsets, int **recursion_offsets);

		int DC_get_blocks_level(int m, int **row_offsets, int **recursion_offsets);

		int* DC_get_RowPerm();

		int* DC_get_RowRev();
							
		DC(int globalNbRow, int partSize, int recursion, int level = 0) 
		{
			treeRoot = new DCNode();
			RowPerm = new int [globalNbRow];
			RowRev = new int [globalNbRow];
			PARTSIZE = partSize;
			LEVEL = level;
			Recursion = recursion;
		}
		
		~DC()
		{
			delete[] RowPerm;
			delete[] RowRev;
			
			delete treeRoot;
		}
};

void DC_create_permutation (int *perm, int *part, int size);

// Permute "tab" 2D array of int using "perm"
void DC_permute (int **tab, int *ntab, int *val, int *rev, int *perm, int nbRow, int offset);
#endif
