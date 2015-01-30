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
	static bool isInside(const std::vector<T> & data, T elem)  { return find(data, elem) == data.size(); }


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

        int nEntry = V.size();
        if (nEntry == 0)  { tmp_stream << "Null"; }
        for (int i = 0; i < nEntry; i++) {
        	tmp_stream << V[i];
        	if (i != nEntry - 1) { tmp_stream << " "; }
        }
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
