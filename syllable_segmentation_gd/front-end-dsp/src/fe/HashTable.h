/*
    This is a header file for HashTable.c

    Copyright (C) 2007-2016 Speech and Music Technology Lab,
    Indian Institute of Technology Madras
    
    Contributed by Chaitanya

    This file is part of KWS.

    KWS is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    KWS is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with KWS.  If not, see <http://www.gnu.org/licenses/>. 
*/


#ifndef HASHTABLE1_H
#define HASHTABLE1_H
#include<stdlib.h>

#define DEFAULT_HASH_TABLE_SIZE 200

typedef struct HashElement
{
	char* data;
	void* fnptr;
	int numFVectors;
	struct HashElement* next;
}HashElement;

typedef struct hashTable
{
        struct HashElement** element;
        int size;
}hashTable;

HashElement* GetHashElement(char* data,void* fnptr,int numFVectors);
void InsertFeaturesInfo(hashTable* ht,int numCepstrum,int numFilters,int fftSize,int windowSize,int resGdWindowSize, int numPthCoefficients);
hashTable* InitHashTable(int size);
void InsertHE(hashTable* ht,char* data,void* fnptr,int numFVectors);
int GetHashValue(int table_size,char* str);
HashElement* SearchHE(struct hashTable *ht,char* data);
void UnloadHashTable(struct hashTable* ht);
#endif
