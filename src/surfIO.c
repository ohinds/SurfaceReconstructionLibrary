/*****************************************************************************
 * surfIO.c is the source file for surface input/output functions for
 * libsr
 * Oliver Hinds <oph@bu.edu> 2005-06-22
 *
 *
 *
 *****************************************************************************/

#include"surfIO.h"

#define SURF_IO_VERSION_C "$Id: surfIO.c,v 1.21 2007/05/22 19:18:11 oph Exp $"

/**
 * reads a surface as a given format
 */
surface *readSurface(char *filename, enum SURFACE_FORMAT format) {
  switch(format) {
    case OFF:
      return readOFF(filename);
      break;
    case OBJ:
      return readOBJ(filename);
      break;
    case MGHSURF:
      return readMGHsurf(filename);
      break;
    default:
      fprintf(stderr,"readSurface(): unsupported format\n");
      break;
  }

  return NULL;
}

/**
 * writes a surface as a given format
 */
int writeSurface(surface *surf, char *filename, enum SURFACE_FORMAT format) {
  switch(format) {
    case OFF:
      return writeOFF(surf,filename);
      break;
    case OBJ:
      return writeOBJ(surf,filename);
      break;
    case MGHSURF:
      return writeMGHsurf(surf,filename);
      break;
    default:
      fprintf(stderr,"writeSurface(): unsupported format\n");
      break;
  }

  return SR_FAILURE;
}

/**
 * reads a surface from an off file
 * returns a surface structure or null if reading flibvps
 */
surface *readOFF(char* filename) {
  int numVerts, numFaces, i, f1, f2, f3;
  double vx, vy, vz;
  FILE *fp;
  surface *surf;

  /* open the file */
  fp = fopen(filename,"r");
  if(fp == NULL) {
    fprintf(stderr,"error: couldn't open off file %s for reading.\n",filename);
    return NULL;
  }

  /* read the header */
  skipComments(fp,'#');
  fscanf(fp, "OFF\n%d %d %d\n", &numVerts, &numFaces, &i);

  if(SR_VERBOSE) {
    fprintf(stdout,"reading off file %s with %d vertices and %d faces...",
            filename,numVerts,numFaces);
  }

  /* create a surface */
  surf = createSurface(numVerts,numFaces);

  /* read the vertices */
  for(i = 0; i < numVerts; i++) {
    fscanf(fp,"%lf",&vx);
    fscanf(fp,"%lf",&vy);
    fscanf(fp,"%lf",&vz);

    /* add the vertex to the surface */
    addVertexCoord(surf,vx,-vy,vz,-1,-1);
  }

  /* read the faces */
  for(i = 0; i < numFaces; i++) {
    fscanf(fp, "%d", &f1);
    fscanf(fp, "%d", &f1);
    fscanf(fp, "%d", &f2);
    fscanf(fp, "%d", &f3);

    /* add the face */
    addFaceInd(surf,f1,f2,f3);
  }
  fclose(fp);

  if(SR_VERBOSE) {
    fprintf(stdout,"done\n");
  }

  return surf;
}

/**
 * writes a surface in the off file formate
 * returns SUCCESS or FAILURE
 */
int writeOFF(surface *surf, char *filename) {
  int i;
  FILE *fp = fopen(filename,"w+");
  if(fp == NULL) {
    fprintf(stderr,"error: surfIO.c couldnt open %s for writing an off file\n",filename);
    return SR_FAILURE;
  }

  /* validate surface */
  if(surf == NULL) {
    fprintf(stderr,"error: can't write off file of NULL surf struct\n");
    return SR_FAILURE;
  }

  if(SR_VERBOSE) {
    fprintf(stdout,"writing off file %s with %d vertices and %d faces...",
            filename,surf->numVertices,surf->numFaces);
  }

  /* write the header */
  fprintf(fp,"OFF\n%d %d 0\n",surf->numVertices,surf->numFaces);

  /* write all the vertices */
  for(i = 0; i < surf->numVertices; i++) {
    fprintf(fp,"%f %f %f\n",
            surf->vertices[i][SR_X],
            surf->vertices[i][SR_Y],
            surf->vertices[i][SR_Z]);
  }

  /* write all the faces */
  for(i = 0; i < surf->numFaces; i++) {
    fprintf(fp,"3 %d %d %d\n",
            surf->faces[i][0],
            surf->faces[i][1],
            surf->faces[i][2]);
  }

  fclose(fp);

  if(SR_VERBOSE) {
    fprintf(stdout,"done\n");
  }

  return SR_SUCCESS;
}

/**
 * reads a surface from an obj file
 * returns a surface structure or null if reading flibvps
 */
surface *readOBJ(char* filename) {
  surface *surf;
  int nv, nf, i;
  double vx,vy,vz;
  int f1,f2,f3;
  char str[10*SR_MAX_STR_LEN];

  /* open the file */
  FILE *fp = fopen(filename,"r");

  /* validate */
  if(fp == NULL) {
    fprintf(stderr,"error: surfIO.c couldn't open the obj file %s for reading\n", filename);
    return NULL;
  }

  skipComments(fp,'#');

  /* take a pass through counting the vertices and faces */
  readLine(fp,str);
  for(nv = 0; str[0] == 'v'; nv++) {
    skipComments(fp,'#');
    readLine(fp,str);
  }

  for(nf = 0; str[0] == 'f'; nf++) {
    skipComments(fp,'#');
    readLine(fp,str);
  }

  /* create a surface */
  surf = createSurface(nv,nf);

  /* validate */
  if(surf == NULL) {
    fprintf(stderr,"error: surfIO.c couldnt allocate memory for a surface\n");
  }

  if(SR_VERBOSE) {
    fprintf(stdout,"reading obj file %s with %d vertices and %d faces...",
            filename,nv,nf);
  }

  /* reset the file pointer at the begining and read the verts and faces */
  fseek(fp,0,0);
  skipComments(fp,'#');

  /* read the vertices */
  for(i = 0; i < nv; i++) {
    fscanf(fp,"v %lf %lf %lf\n",&vx,&vy,&vz);
    addVertexCoord(surf,vx,-vy,vz,-1,-1);
    skipComments(fp,'#');
  }

  /* read the faces */
  for(i = 0; i < nf; i++) {
    fscanf(fp,"f %d %d %d\n", &f1, &f2, &f3);
    addFaceInd(surf,f1-1,f2-1,f3-1);
    skipComments(fp,'#');
  }

  fclose(fp);

  if(SR_VERBOSE) {
    fprintf(stdout,"done\n");
  }

  return surf;
}

/**
 * write an obj file
 */
int writeOBJ(surface *surf, char *filename) {
  int i;
  FILE *fp = fopen(filename,"w+");
  if(fp == NULL) {
    fprintf(stderr,"error: surfIO.c couldnt open %s for writing an obj file\n",filename);
    return SR_FAILURE;
  }

  /* validate surface */
  if(surf == NULL) {
    fprintf(stderr,"error: can't write obj file of NULL surf struct\n");
    return SR_FAILURE;
  }

  /* write the header */
  if(SR_VERBOSE) {
    fprintf(stderr,"writing obj with %d vertices and %d faces\n",
            surf->numVertices, surf->numFaces);
  }

  /* write the header */
  fprintf(fp,"# generated by surfRecon\n");

  /* write all the vertices */
  for(i = 0; i < surf->numVertices; i++) {
    fprintf(fp,"v %f %f %f\n",
            surf->vertices[i][SR_X],
            surf->vertices[i][SR_Y],
            surf->vertices[i][SR_Z]);
  }

  /* write all the faces */
  for(i = 0; i < surf->numFaces; i++) {
    fprintf(fp,"f %d %d %d\n",
            surf->faces[i][0]+1,
            surf->faces[i][1]+1,
            surf->faces[i][2]+1);
  }

  fclose(fp);

  if(SR_VERBOSE) {
    fprintf(stdout,"done\n");
  }

  return SR_SUCCESS;
}

/**
 * reads a surface in a freesurfer triangle surface binary format
 * returns a surface structure or null if reading fails
 *
 * got file format from
 * http://wideman-one.com/gw/brain/fs/surfacefileformats.htm
 *
 * used a couple tricks from mukund's readMGHsurf.m
 */
#define TRIANGLE_FILE_MAGIC_NUMBER 16777214
#define QUAD_FILE_MAGIC_NUMBER 16777215
#define NEW_QUAD_FILE_MAGIC_NUMBER 16777213
surface *readMGHsurf(char* filename) {
  unsigned int magicNum, numVerts, numFaces, i, f[4];
  unsigned short s[3];
  char trash[SR_MAX_STR_LEN];
  FILE *fp;
  float v[3];
  surface *surf = NULL;

  /* open the file */
  fp = fopen(filename,"r");
  if(fp == NULL) {
    fprintf(stderr,"error: couldn't open fs surface file %s for reading.\n",filename);
    return NULL;
  }

  /* read the header */
  readBigEndianInt24(fp,1,&magicNum);

  /* test the magic number */
  if(magicNum != TRIANGLE_FILE_MAGIC_NUMBER
     && magicNum != QUAD_FILE_MAGIC_NUMBER
     && magicNum != NEW_QUAD_FILE_MAGIC_NUMBER
     ) {
    fprintf(stderr,"error: magic number (%u) doesn't match any MGH surface formats for %s.\n", magicNum, filename);
    return NULL;
  }

  if(magicNum == TRIANGLE_FILE_MAGIC_NUMBER) {
    /* two '\n's */
    readLine(fp,trash);
    readLine(fp,trash);

    /* read the sizes */
    readBigEndianInt32(fp,1,&numVerts);
    readBigEndianInt32(fp,1,&numFaces);

    surf = createSurface(numVerts,numFaces);

    /* read the vertices */
    for(i = 0; i < numVerts; i++) {
      readBigEndianFloat32(fp,3,v);
      addVertexCoord(surf,v[0],v[1],v[2],-1,-1);
    }

    /* read the faces */
    for(i = 0; i < numFaces; i++) {
      readBigEndianInt32(fp,3,f);
      addFaceInd(surf,f[0],f[1],f[2]);
    }
  }
  else if(magicNum == QUAD_FILE_MAGIC_NUMBER) {
    /* read the sizes */
    readBigEndianInt24(fp,1,&numVerts);
    readBigEndianInt24(fp,1,&numFaces);

    surf = createSurface(numVerts,2*numFaces);

    /* read the vertices */
    for(i = 0; i < numVerts; i++) {
      readBigEndianShort16(fp,3,s);
      addVertexCoord(surf,s[0]/100.0f,v[1]/100.0f,v[2]/100.0f,-1,-1);
    }

    /* read the faces */
    for(i = 0; i < numFaces; i++) {
      readBigEndianInt24(fp,4,f);
      addFaceInd(surf,f[0],f[1],f[2]);
      addFaceInd(surf,f[0],f[2],f[3]);
    }
  }
  else if(magicNum == NEW_QUAD_FILE_MAGIC_NUMBER) {
    /* read the sizes */
    readBigEndianInt24(fp,1,&numVerts);
    readBigEndianInt24(fp,1,&numFaces);

    surf = createSurface(numVerts,2*numFaces);

    /* read the vertices */
    for(i = 0; i < numVerts; i++) {
      readBigEndianFloat32(fp,3,v);
      addVertexCoord(surf,v[0],v[1],v[2],-1,-1);
    }

    /* read the faces */
    for(i = 0; i < numFaces; i++) {
      readBigEndianInt24(fp,4,f);
      addFaceInd(surf,f[0],f[1],f[2]);
      addFaceInd(surf,f[0],f[2],f[3]);
    }
  }
  else {
    // wont get here
    return NULL;
  }

  fclose(fp);

  fprintf(stderr,"read MGH surf file with %d vertices.\n",numVerts);

  return surf;
}

/**
 * write a surface in a freesurfer triangle surface binary format
 */
int writeMGHsurf(surface *surf, char *filename) {
  int i;
  unsigned int magicNum = TRIANGLE_FILE_MAGIC_NUMBER;
  float f[3];
  char *tm;
  FILE *fp = fopen(filename,"w+");
  if(fp == NULL) {
    fprintf(stderr,"error: surfIO.c couldnt open %s for writing a fs surf\n",filename);
    return SR_FAILURE;
  }

  /* validate surface */
  if(surf == NULL) {
    fprintf(stderr,"error: can't write fs surf of NULL surf struct\n");
    return SR_FAILURE;
  }

  fprintf(stderr,"writing freesurfer binary triangle surface with %d vertices and %d faces\n",
          surf->numVertices, surf->numFaces);

  /* write the header */
  writeBigEndianInt24(fp,1,&magicNum);

  tm = getTimeString();
  fprintf(fp,"created by libsr on %s\n\n",tm);
  free(tm);

  writeBigEndianInt32(fp,1,&surf->numVertices);
  writeBigEndianInt32(fp,1,&surf->numFaces);

  /* write the vertices */
  for(i = 0; i < surf->numVertices; i++) {
    f[0] = (float) surf->vertices[i][0];
    f[1] = (float) surf->vertices[i][1];
    f[2] = (float) surf->vertices[i][2];
    writeBigEndianFloat32(fp,3,f);
  }

  /* write the faces */
  for(i = 0; i < surf->numFaces; i++) {
    writeBigEndianInt32(fp,3,surf->faces[i]);
  }

  fclose(fp);
  return SR_SUCCESS;
}

/**
 * read a label file
 */
int readLabelFile(surface *surf, char *filename) {
  int i;
  FILE *fp = fopen(filename,"r");
  if(fp == NULL) {
    fprintf(stderr,"error: surfIO.c couldnt open %s for reading a label file\n",filename);
    return SR_FAILURE;
  }

  /* validate surface */
  if(surf == NULL || surf->vertexLabels == NULL) {
    fprintf(stderr,"error: can't read label file into NULL surf struct\n");
    return SR_FAILURE;
  }

  /* zero the label */
  for(i = 0; i < surf->maxVertices; i++) {
    surf->vertexLabels[i] = -1;
  }

  skipComments(fp,'#');

  /* read all the vertex labels */
  for(i = 0; !feof(fp) && i < surf->numVertices; i++) {
    fscanf(fp,"%lf",&surf->vertexLabels[i]);
  }

  fclose(fp);

  return SR_SUCCESS;
}

/**
 * write a label file
 */
int writeLabelFile(surface *surf, char *filename) {
  int i;
  FILE *fp = fopen(filename,"w+");
  if(fp == NULL) {
    fprintf(stderr,"error: surfIO.c couldnt open %s for writing a label file\n",filename);
    return SR_FAILURE;
  }

  /* validate surface */
  if(surf == NULL || surf->vertexLabels == NULL) {
    fprintf(stderr,"error: can't write label file of NULL surf struct\n");
    return SR_FAILURE;
  }

  /* write the header */
  fprintf(fp,"# generated by surfRecon\n");

  /* write all the vertex labels */
  for(i = 0; i < surf->numVertices; i++) {
    fprintf(fp,"%lf\n",surf->vertexLabels[i]);
  }

  fclose(fp);

  return SR_SUCCESS;
}

/**
 * read a label file
 */
int readCurvatureFile(surface *surf, char *filename) {
  int i;
  FILE *fp = fopen(filename,"r");

  if(fp == NULL) {
    fprintf(stderr,"error: surfIO.c couldnt open %s for reading a curvature file\n",filename);
    return SR_FAILURE;
  }

  /* validate surface */
  if(surf == NULL) {
    fprintf(stderr,"error: can't read curvature file into NULL surf struct\n");
    return SR_FAILURE;
  }

  /* create curvature, if needed */
  if(surf->curvature == NULL) {
    surf->curvature = (int*) malloc(sizeof(int)*max(surf->numFaces,surf->numVertices));
  }

  skipComments(fp,'#');

  /* read all the curvature data */
  for(i = 0; !feof(fp) && i < max(surf->numVertices,surf->numFaces); i++) {
    fscanf(fp,"%d",&surf->curvature[i]);
  }

  /* set the mode */
  if(abs(i-surf->numFaces) < abs(i-surf->numVertices)) {
    surf->curvatureMode = FACE_CURVATURE;
  }
  else {
    surf->curvatureMode = VERTEX_CURVATURE;
  }

  fclose(fp);

  return SR_SUCCESS;
}

/**
 * write a label file
 */
int writeCurvatureFile(surface *surf, char *filename) {
  int i;
  FILE *fp = fopen(filename,"w+");
  if(fp == NULL) {
    fprintf(stderr,"error: surfIO.c couldnt open %s for writing a label file\n",filename);
    return SR_FAILURE;
  }

  /* validate surface */
  if(surf == NULL || surf->curvature == NULL || surf->curvatureMode == NO_CURVATURE) {
    fprintf(stderr,"error: can't write curvature file of NULL surf struct\n");
    return SR_FAILURE;
  }

  /* write the header */
  fprintf(fp,"# generated by surfRecon\n");

  /* write all the curvature labels */
  for(i = 0; i < (surf->curvatureMode == VERTEX_CURVATURE ? surf->numVertices : surf->numFaces); i++) {
    fprintf(fp,"%d\n",surf->curvature[i]);
  }

  fclose(fp);

  return SR_SUCCESS;
}

/** freesurfer label and curvature formats **/

/**
 * read a freesurfer label file
 */
int readMGHLabelFile(surface *surf, char *filename) {
  int i, v, numLabels = 0;
  float tr;
  double lab;
  char str[SR_MAX_STR_LEN] = "";
  FILE *fp = fopen(filename,"r");
  if(fp == NULL) {
    fprintf(stderr,"error: surfIO.c couldnt open %s for reading an MGH label file\n",filename);
    return SR_FAILURE;
  }

  /* validate surface */
  if(surf == NULL || surf->vertexLabels == NULL) {
    fprintf(stderr,"error: can't read label file into NULL surf struct\n");
    return SR_FAILURE;
  }

  /* zero the label */
  for(i = 0; i < surf->maxVertices; i++) {
    surf->vertexLabels[i] = -1;
  }

  skipComments(fp,'#');

  /* read the number of labeled vertices */
  fscanf(fp,"%d",&numLabels);

  /* read all the vertex labels */
  for(i = 0; !feof(fp) && i < numLabels; i++) {
    fscanf(fp,"%d %f %f %f %lf",&v,&tr,&tr,&tr,&lab);
    readLine(fp,str);


    if(v < surf->numVertices) {
      surf->vertexLabels[v] = lab;
    }
    else {
      fprintf(stderr,"vertex index %d out of bounds. num vertices is %d\n",
              v, surf->numVertices);
    }
  }

  fclose(fp);

  return SR_SUCCESS;
}

/**
 * write a freesurfer label file
 */
int writeMGHLabelFile(surface *surf, char *filename) {
  return notSupported();
}

/**
 * read a freesurfer curvature file
 */
#define CURVATURE_FILE_MAGIC_NUMBER 16777215
int readMGHCurvatureFile(surface *surf, char *filename) {
  unsigned int i,numVerts,numFaces,magicNum;
  float c;
  //float v[3];
  FILE *fp = fopen(filename,"r");

  if(fp == NULL) {
    fprintf(stderr,"error: surfIO.c couldnt open %s for reading an MGH curvature file\n",filename);
    return SR_FAILURE;
  }

  /* validate surface */
  if(surf == NULL) {
    fprintf(stderr,"error: can't read curvature file into NULL surf struct\n");
    fclose(fp);
    return SR_FAILURE;
  }

  /* create curvature, if needed */
  if(surf->curvature == NULL) {
    surf->curvature = (int*) malloc(sizeof(int)*surf->numVertices);
  }

  /* read the header */
  readBigEndianInt24(fp,1,&magicNum);

  /* test the magic number */
  if(magicNum != CURVATURE_FILE_MAGIC_NUMBER) {
    fprintf(stderr,"error: magic number (%u) doesn't match (%u) for %s.\n", magicNum, CURVATURE_FILE_MAGIC_NUMBER, filename);
    return SR_FAILURE;
  }

  /* read surf sizes and see if they match */
  readBigEndianInt32(fp,1,&numVerts);
  readBigEndianInt32(fp,1,&numFaces);
  readBigEndianInt32(fp,1,&i);

  if(numVerts != surf->numVertices || numFaces != surf->numFaces) {
    fprintf(stderr,"error: surface and  curvature file don't match\n");
    fclose(fp);
    return SR_FAILURE;
  }

  /* read all the curvature data */
  for(i = 0; !feof(fp) && i < surf->numVertices; i++) {
    //readBigEndianFloat32(fp,3,v);
    readBigEndianFloat32(fp,1,&c);
    surf->curvature[i] = (c < 0 ? -1 : 1);
  }

  /* set the mode */
  surf->curvatureMode = VERTEX_CURVATURE;

  fclose(fp);

  return SR_SUCCESS;
}

/**
 * write a freesurfer curvature file
 */
int writeMGHCurvatureFile(surface *surf, char *filename) {
  return notSupported();
}

/**
 * read a set of slice contours from a file
 */
list *readSliceContours(char *filename) {
  list *slices, *slice, *nextSlice;
  contour *cont;
  listNode *sliceNode, *contNode, *adjNode, *lnAdj;
  intNode *in;
  vertex *v;
  char str[10*SR_MAX_STR_LEN];
  double sliceZ;
  int adjIndex;

  /* open the file */
  FILE *fp = fopen(filename,"r");

  /* validate */
  if(fp == NULL) {
    fprintf(stderr,"error: surfIO.c couldn't open the slice contour file %s for reading\n", filename);
    return NULL;
  }

  skipComments(fp,'#');

  /* create the slices */
  slices = newList(LIST);

  /* read the first string */
  fscanf(fp, "slice ");

  /* read the file */
  while(!feof(fp)) {
    /* read the z coordinate */
    fscanf(fp, "%lf", &sliceZ);

    slice = newList(LIST);
    enqueue(slices,slice);

    /* read the next string and test for new slice */
    fscanf(fp, "%s ", str);
    while(!strcmp(str,"contour")) { /* while we are still reading contours */
      cont = createContour();
      enqueue(slice,cont);

      /* read closure state */
      fscanf(fp, "%s", str);
      if(!strcmp(str,"closed")) {
        cont->closed = CLOSED;
      }

      /* read the first string on the next line, test for "adjacent" */
      fscanf(fp, "%s", str);

      while(strcmp(str,"adjacent")) {
        /* create and read this vertex */
        v = createVertex();
        v->x = atof(str);
        fscanf(fp, "%lf", &v->y);
        v->z = sliceZ;

        enqueue(cont->vertices,v);

        fscanf(fp, "%s", str);
      }

      /* read the adjacent contours as ints, convert to pointers later */
      fscanf(fp, "%d", &adjIndex);
      while(adjIndex != -1) {
        /* add this index to the list */
        in = (intNode*) malloc(sizeof(intNode));
        in->val = adjIndex;

        enqueue(cont->adjacentContours,in);

        fscanf(fp, "%d", &adjIndex);
      }

      fscanf(fp, "%s ", str);
    }
  }

  /* resolve the links from the ajacent contour indices */
  for(sliceNode = getListNode(slices,0);
      listSize(slices) > 1 &&((listNode*)sliceNode)->next;
      sliceNode = (listNode*) sliceNode->next) {

    /* get links to this slice list and the next, for searching */
    slice = (list*) sliceNode->data;
    nextSlice = (list*) ((listNode*) ((listNode*)sliceNode)->next)->data;

    /* iterate over the contours in this slice */
    for(contNode = getListNode(slice,0); contNode;
        contNode = (listNode*) contNode->next) {
      cont = (contour*) contNode->data;

      /* iterate over the adjacent contour indices */
      for(adjNode = getListNode(cont->adjacentContours,0); adjNode;
          adjNode = (listNode*) adjNode->next) {

        /* get the index, replace it with a link, or NULL */
        in = (intNode*) adjNode->data;
        lnAdj = getListNode(nextSlice,in->val);
        free(in);

        if(lnAdj == NULL) {
          markForDeletion(adjNode);
        }
        else {
          setListNodeData(adjNode, (void*) lnAdj->data);
        }
      }

      deleteMarkedNodes(cont->adjacentContours);
    }
  }

  fclose(fp);

  return slices;
}

/**
 * write a set of slice contours to a file
 */
int writeSliceContours(list *slices, char *filename) {
  listNode *sliceNode, *contNode, *vertNode, *adjNode;
  contour *cont;
  vertex *v;
  int i;

  /* open the file */
  FILE *fp = fopen(filename,"w");


  /* validate */
  if(fp == NULL) {
    fprintf(stderr,"error: surfIO.c couldn't open the slice contour file %s for writing\n", filename);
    return SR_FAILURE;
  }

  //fprintf(fp,"# slice contour file created by libsr\n# %s\n\n", SURF_IO_VERSION_C);

  /* write each slice */
  for(sliceNode = getListNode(slices,0), i = 0; sliceNode;
      sliceNode = (listNode*) sliceNode->next, i++) {
    if(sliceEmpty((list*) sliceNode->data)) continue;

    fprintf(fp,"slice %lf\n",
            ((vertex*) getListNode(((contour*)
                                    getListNode((list*) sliceNode->data,
                                                0)->data)->vertices,
                                   0)->data)->z);

    /* iterate over the contours in this slice */
    for(contNode = getListNode((list*) sliceNode->data,0); contNode;
        contNode = (listNode*) contNode->next) {
      cont = (contour*) contNode->data;

      /* print the info for this contour */
      fprintf(fp, "contour ");
      switch(cont->closed) {
        case OPEN:
          fprintf(fp,"open");
          break;
        case CLOSED:
          fprintf(fp,"closed");
          break;
      }
      fprintf(fp,"\n");

      /* iterate over the vertices, writing each */
      for(vertNode = getListNode(cont->vertices,0); vertNode;
          vertNode = (listNode*) vertNode->next) {
        v = (vertex*) vertNode->data;

        fprintf(fp,"%lf %lf\n", v->x, v->y);
      }


      fprintf(fp,"adjacent\n");
      /* iterate over the adjacent contours, writing the indices */
      for(adjNode = getListNode(cont->adjacentContours,0);
          i+1 < listSize(slices)&& adjNode;
          adjNode = (listNode*) adjNode->next) {
        fprintf(fp, "%d ", findInListI((list*) getListNode(slices,i+1)->data,
                                       adjNode->data));
      }
      fprintf(fp,"-1\n");

    }
  }

  fclose(fp);

  return SR_SUCCESS;
}

/**
 * write a set of slice contour labels to a file
 */
int writeSliceContourLabels(list *slices, char *filename) {
  listNode *sliceNode, *contNode, *vertNode;
  contour *cont;
  vertex *v;
  int i;

  /* open the file */
  FILE *fp = fopen(filename,"w");

  /* validate */
  if(fp == NULL) {
    fprintf(stderr,"error: surfIO.c couldn't open the slice contour file %s for writing\n", filename);
    return SR_FAILURE;
  }

  //fprintf(fp,"# slice contour file created by libsr\n# %s\n\n", SURF_IO_VERSION_C);

  /* write each slice */
  for(sliceNode = getListNode(slices,0), i = 0; sliceNode;
      sliceNode = (listNode*) sliceNode->next, i++) {
    /* iterate over the contours in this slice */
    for(contNode = getListNode((list*) sliceNode->data,0); contNode;
        contNode = (listNode*) contNode->next) {
      cont = (contour*) contNode->data;

      /* iterate over the vertices, writing each label */
      for(vertNode = getListNode(cont->vertices,0); vertNode;
          vertNode = (listNode*) vertNode->next) {
        v = (vertex*) vertNode->data;

        fprintf(fp,"%d\n", v->label);
      }
    }
  }

  fclose(fp);

  return SR_SUCCESS;
}

/**
 * read a set of slice contour labels from a file
 */
int readSliceContourLabels(list *slices, char *filename) {
  listNode *sliceNode, *contNode, *vertNode;
  contour *cont;
  vertex *v;
  int i;

  /* open the file */
  FILE *fp = fopen(filename,"r");

  /* validate */
  if(fp == NULL) {
    fprintf(stderr,"error: surfIO.c couldn't open the slice contour file %s for reading\n", filename);
    return SR_FAILURE;
  }

  //fprintf(fp,"# slice contour file created by libsr\n# %s\n\n", SURF_IO_VERSION_C);

  /* read each label */
  for(sliceNode = getListNode(slices,0), i = 0; !feof(fp) && sliceNode;
      sliceNode = (listNode*) sliceNode->next, i++) {
    /* iterate over the contours in this slice */
    for(contNode = getListNode((list*) sliceNode->data,0); !feof(fp) && contNode;
        contNode = (listNode*) contNode->next) {
      cont = (contour*) contNode->data;

      /* iterate over the vertices, writing each label */
      for(vertNode = getListNode(cont->vertices,0); !feof(fp) && vertNode;
          vertNode = (listNode*) vertNode->next) {
        v = (vertex*) vertNode->data;

        fscanf(fp,"%d\n", &v->label);
      }
    }
  }

  fclose(fp);

  return SR_SUCCESS;
}

/**
 * delete a list of slice contours
 */
void deleteSliceContours(list *slices) {
  listNode *i,*j;
  int k,l;
  list *slice;
  contour *cont;
  for(i = getListNode(slices,0), k = 0; i; k++, i = (listNode*) i->next) {
    slice = (list*) i->data;
    for(j = getListNode(slice,0), l = 0; j; l++, j = (listNode*) j->next) {
      cont = (contour*) j->data;
      deleteContour(cont);
    }
    freeList(slice);
  }
  freeList(slices);
}

/**
 * reads a nuages file into a set of slices
 */
list *readNuagesSlices(char *filename) {
  FILE *fp;
  list *slices, *slice;
  contour *cont;
  vertex *vert;
  double curZ;
  int numSlices, numVerts;
  char str[SR_MAX_STR_LEN] = "";

  // open the file
  if(NULL == (fp = fopen(filename,"r"))) {
    return NULL;
  }

  // allocate
  slices = newList(LIST);
  slice = newList(LIST);
  enqueue(slices,slice);
  cont = createContour();
  cont->closed = CLOSED;
  enqueue(slice,cont);
  fscanf(fp, "S %d\nv %d z %lf\n{", &numSlices, &numVerts, &curZ);

  // read all contours
  while(!feof(fp)) {
    fscanf(fp, "%s", str);
    if(!strcmp(str,"}")) {
      fscanf(fp, "%s", str);
      if(feof(fp)) break; // if we read past the end

      if(!strcmp(str,"v")) { // new slice
        fscanf(fp, "%d z %lf\n{",&numVerts,&curZ);
        slice = newList(LIST);
        enqueue(slices,slice);
      }
      cont = createContour();
      cont->closed = CLOSED;
      enqueue(slice,cont);
    }
    else {
      // create a new vertex
      vert = createVertex();
      vert->x = atof(str);
      fscanf(fp,"%lf\n",&vert->y);
      vert->z = curZ;
      enqueue(cont->vertices,vert);
    }
  }

  return slices;
}

/**
 * looks for a valid surface format name in a string,
 * returns both the format
 */
int getSurfFormatFromString(char *str, enum SURFACE_FORMAT *format) {
  /* validate */
  if(str == NULL || str[0] == '\0') {
    *format = INVALID_SURFACE_FORMAT;
    return SR_FAILURE;
  }

  /* look for a valid type */
  if(!strcmp(str,"off") || !strcmp(str,"OFF")) {
    *format = OFF;
  }
  else if(!strcmp(str,"obj") || !strcmp(str,"OBJ")) {
    *format = OBJ;
  }
  else if(!strcmp(str,"mgh") || !strcmp(str,"MGH") ||
          !strcmp(str,"mghsurf") || !strcmp(str,"MGHsurf") ||
          !strcmp(str,"mghSurf") || !strcmp(str,"MGHSurf") ||
          !strcmp(str,"MGHSURF") || !strcmp(str,"fs") || !strcmp(str,"FS")) {
    *format = MGHSURF;
  }
  else {
    *format = INVALID_SURFACE_FORMAT;
    return SR_FAILURE;
  }

  return SR_SUCCESS;
}

/********************************************************************
 * $Source: /home/cvs/PROJECTS/SurfaceReconstructionLibrary/src/surfIO.c,v $
 * Local Variables:
 * mode: C
 * fill-column: 76
 * comment-column: 0
 * End:
 ********************************************************************/
