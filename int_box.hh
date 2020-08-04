#ifndef INT_BOX_HH
#define INT_BOX_HH

struct int_box {
    /** The lower x bound (inclusive). */
    int li;
    /** The upper x bound (exclusive). */
    int ui;
    /** The lower y bound (inclusive). */
    int lj;
    /** The upper y bound (exclusive). */
    int uj;
    /** Creates an integer box type, initializing the four index bounds.
     * \param[in] (li_,ui_) the x index bounds.
     * \param[in] (lj_,uj_) the y index bounds. */
    int_box(int li_,int ui_,int lj_,int uj_) :
        li(li_), ui(ui_), lj(lj_), uj(uj_) {}
    /** Creates an integer box type, initializing the four index bounds
     * from a pointer to memory.
     * \param[in] pp the pointer to read from. */
    int_box(int *pp) : li(*pp), ui(pp[1]), lj(pp[2]), uj(pp[3]) {}
    /** Extends the bounds so that the include the given grid point.
     * \param[in] (i,j) the point to consider. */
    inline void bound(int i,int j) {
        if(i<li) li=i;
        if(i>=ui) ui=i+1;
        if(j<lj) lj=j;
        if(j>=uj) uj=j+1;
    }
    /** Extends the data structure by an equal padding amount in each
     * direction.
     * \param[in] pad the number of grid points to pad by. */
    inline void extend(int pad) {
        li-=pad;ui+=pad;
        lj-=pad;uj+=pad;
    }
    /** Trims the box size to lie within the given bounds.
     * \param[in] (li_,ui_) the x bounds.
     * \param[in] (lj_,uj_) the y bounds. */
    inline void trim(int li_,int ui_,int lj_,int uj_) {
        if(li<li_) li=li_;
        if(ui>ui_) ui=ui_;
        if(lj<lj_) lj=lj_;
        if(uj>uj_) uj=uj_;
    }
    /** Computes the width of the box in the x direction.
     * \return The width. */
    inline int isize() {return ui-li;}
    /** Computes the height of the box in the y direction.
     * \return The height. */
    inline int jsize() {return uj-lj;}
    /** Computes the total number of grid points covered by the box.
     * \return The number of grid points. */
    inline int size() {return isize()*jsize();}
};

#endif
