#ifndef DUNE_VECTORHELPER_HH
#define DUNE_VECTORHELPER_HH

#include <limits>
#include <map>
#include <string>
#include <utility>
#include <vector>
#include <iostream>
#include <assert.h>

#include <fstream>


namespace Dune {


class VectorHelper
{

public:

	VectorHelper()  { }


	// Finds an element inside sorted non-repeating vector, returns index of the element if found
	// If not found, returns size of the vector
	//
	// [TODO] Implement binary search, but atm for small arrays direct search is fine
	template<typename T>
	static int find(const std::vector<T> & data, T elem)
	{
		for (int i = 0; i < data.size(); i++)  { if (data[i] == elem)  { return i; }}
		return data.size();
	}


	// Checks if an element is inside a sorted non-repeating vector
	template<typename T>
	static bool isInside(const std::vector<T> & data, T elem)  { return find(data, elem) < data.size(); }


	// Sorts array and removes all repeating entries. No overhead
	template<typename T>
	static void compactify(std::vector<T> & data)
	{
		if (data.size() > 0)
		{
			std::sort(data.begin(), data.end());

			// Delete all repeating data
			int j = 1;
			for (int i = 1; i < data.size(); i++)
			{
				if (data[i-1] != data[i]) { data[j++] = data[i]; }
			}

			data.resize(j);
		}
	}


    // Converts an arbitrary vector into string by sticking all of the entries together
    // Whitespace separated
	template<typename T>
    static std::string vector2string(const T & V)
    {
        std::stringstream tmp_stream;
        tmp_stream << "(";

        int nEntry = V.size();
        if (nEntry == 0)  { tmp_stream << "Null"; }
        for (int i = 0; i < nEntry; i++) {
        	tmp_stream << V[i];
        	if (i != nEntry - 1) { tmp_stream << ", "; }
        }
        tmp_stream << ")";
        return tmp_stream.str();
    }


    // Converts an arbitrary vector into string by sticking all of the entries together
    // Whitespace separated
	template<typename T>
    static std::string map2string(const T & M)
    {
		typedef typename T::const_iterator  MapIter;

        std::stringstream tmp_stream;
        tmp_stream << "(";

        int nEntry = M.size();
        if (nEntry == 0)  { tmp_stream << "Null"; }

        MapIter iterB = M.begin();
        MapIter iterE = M.end();

        for (MapIter iter = iterB; iter != iterE; iter++) {
        	if (iter != iterB) { tmp_stream << ", "; }
        	tmp_stream << "[" << (*iter).first << ";" << (*iter).second << "]";
        }
        tmp_stream << ")";
        return tmp_stream.str();
    }


    // Takes two sorted arrays with non-repeating entries
    // Returns an array which only has entries found in both input arrays
	template<typename T>
    static std::vector<T> sortedSetIntersection(const std::vector<T> & A, const std::vector<T> & B)
    {
        std::vector<T>  rez;

        int indA = 0;
        int indB = 0;

        while ((indA < A.size()) && (indB < B.size()))
        {
                  if (A[indA] < B[indB]) { indA++; }
             else if (A[indA] > B[indB]) { indB++; }
             else {
                 rez.push_back(A[indA]);
                 indA++;
                 indB++;
             }
        }

        return rez;
    }


    // Takes two sorted arrays with non-repeating entries
    // Returns an array which only has entries found in both input arrays
	template<typename T>
    static std::vector<T> sortedSetUnion(const std::vector<T> & A, const std::vector<T> & B)
    {
        std::vector<T>  rez;

        int indA = 0;
        int indB = 0;

        while ((indA < A.size()) && (indB < B.size()))
        {
                  if (A[indA] < B[indB]) { rez.push_back(A[indA++]); }
             else if (A[indA] > B[indB]) { rez.push_back(B[indB++]); }
             else {
                 rez.push_back(A[indA]);
                 indA++;
                 indB++;
             }
        }

        // If one array runs out before the other, add all remaining entries to the result as they are not duplicate
        for (int i = indA; i < A.size(); i++)  { rez.push_back(A[i]); }
        for (int i = indB; i < B.size(); i++)  { rez.push_back(B[i]); }

        return rez;
    }


    // Takes two sorted arrays with non-repeating entries
    // Returns array A / B (as a set)
	template<typename T>
	static std::vector<T> sortedSetComplement(const std::vector<T> & A, const std::vector<T> & B)
	{
		std::vector<T> rez;
		int indA = 0;
		int indB = 0;

		while ((indA < A.size())&&(indB < B.size()))
		{
			     if (A[indA] >  B[indB])  { indB++;          }
			else if (A[indA] == B[indB])  { indA++;  indB++; }
			else                          { rez.push_back(A[indA++]); }
		}

		for (int i = indA; i < A.size(); i++) { rez.push_back(A[i]); }

		return rez;
	}


	// For two arrays of the same size, no repeating entries, some entries present in both arrays, output an array of indices
	// The output is indexed by the position in the original array A, and points to the position in the original array B
	// If at this position of A there is a shared entity, it should point at corresponding position in B, otherwise at any other position
	// Note: The output should span the entire array B (no repeating indices)
	template<typename T>
	static std::vector<int> partialMatchIndex(const std::vector<T> & A, const std::vector<T> & B)
	{
		typedef std::map<T, int> IndexMap;
		typedef typename IndexMap::iterator  IndexMapIter;

		assert(A.size() == B.size());  // This function is only defined for arrays of the same size

		std::vector<int> unusedB;
		std::vector<int> rez;

		IndexMap AA;
		IndexMap BB;

		// Loop over all elements of A and B, and initialize maps for them
		for (int i = 0; i < A.size(); i++)
		{
			IndexMapIter iterA = AA.find(A[i]);   assert(iterA == AA.end());   // No repeating elements in arrays
			IndexMapIter iterB = BB.find(B[i]);   assert(iterB == BB.end());   // No repeating elements in arrays

			AA.insert(std::pair<T, int>(A[i], i));
			BB.insert(std::pair<T, int>(B[i], i));
		}

		// Loop over all elements of B, and if they are not in A, push their index back to unused set
		// Thus we artificially introduce order to all the entries of B that are not shared with A
		for (IndexMapIter iterB = BB.begin(); iterB != BB.end(); iterB++)
		{
			if (AA.find((*iterB).first) == AA.end())  { unusedB.push_back((*iterB).second); }
		}

		// Loop over all elements in A
		int unusedIndex = 0;
		for (int i = 0; i < A.size(); i++)
		{
			IndexMapIter iterB = BB.find(A[i]);
			if (iterB == BB.end())  { rez.push_back(unusedB[unusedIndex++]); }  //  If they are not in B, use next element from unused set
			else                    { rez.push_back((*iterB).second); }         //  If they are in B, use B'th index
		}

		return rez;
	}


	// Assume sorted array
	// Should use binary search
	/*
	static void remove(std::vector<T> & data, T entry)
	{
		int pos = find(data, entry);
		data.erase(data.begin() + pos);
	}
	*/


};

} // namespace Dune

#endif  // DUNE_VECTORHELPER_HH
