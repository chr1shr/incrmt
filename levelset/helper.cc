#include "helper.hh"

/** This routine prints a fatal error message and exits.
 * \param[in] p The message to print. */
void levelpp_fatal_error(const char *p,const int status) {
    fprintf(stderr,"level++: %s\n",p);
    exit(1);
}

template<class f_type>
void output(const char* filename,f_type *array,int m,int n,double ax,double ay,double dx,double dy) {
    int i,j;
    float *buf=new float[m+1],*fp=buf;
    f_type *arrayp=array,*arraye;
puts(filename);
    // Open the file
    FILE *outf=fopen(filename,"wb");
    if(outf==NULL) levelpp_fatal_error("Can't open output file",1);

    // Save the header line
    *(fp++)=m;
    for(i=0;i<m;i++) *(fp++)=ax+i*dx;
    if(fwrite(buf,sizeof(float),m+1,outf)!=size_t(m+1)) levelpp_fatal_error("Error writing to output file",1);

    // Save the lines of the field
    for(j=0;j<n;j++) {
        fp=(float*) buf;
        *(fp++)=ay+j*dy;
        arraye=arrayp+m;
        while(arrayp<arraye) *(fp++)=*(arrayp++);
        if(fwrite(buf,sizeof(float),m+1,outf)!=size_t(m+1)) levelpp_fatal_error("Error writing to output file",1);
    }

    // Close the file
    fclose(outf);
    delete [] buf;
}

// Explicit instantiation
template void output(const char*,double*,int,int,double,double,double,double);
template void output(const char*,int*,int,int,double,double,double,double);
template void output(const char*,unsigned int*,int,int,double,double,double,double);
