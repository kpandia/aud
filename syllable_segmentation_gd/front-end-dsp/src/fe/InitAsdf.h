/*
    This is a header file for InitAsdf.c

    Copyright (C) 1997-2016 Speech and Music Technology Lab,
    Indian Institute of Technology Madras
    
    Contributed by Hema A Murthy <hema@cse.iitm.ac.in>, Chaitanya

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



#ifndef INIT_ASDF_H
#define INIT_ASDF_H
void PutIAttribute(ASDF *asdf, char *attribute_name, int attribute_value);
void PutFAttribute(ASDF *asdf, char *attribute_name, float attribute_value);
void InitializeASDF (ASDF *asdf);
void InitializeStandardFrontEnd(ASDF *asdf,FILE *fp);
F_VECTOR *GsfRead(ASDF *asdf, int frame_index, char *feature_name);
F_VECTOR *ComputeFeature(ASDF *asdf, int frame_index, char *feature_name);
void GsfOpen(ASDF *asdf, char *filename);
void GsfClose(ASDF *asdf);
int GetIAttribute(ASDF *asdf, char *string);
float GetFAttribute(ASDF *asdf, char *string);
void *GetPtrAttribute(ASDF *asdf, char *string);
int select_channel;
#endif

