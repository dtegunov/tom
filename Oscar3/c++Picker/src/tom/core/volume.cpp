/***********************************************************************//**
 * \file volume.cpp
 * \brief
 * \author  Thomas Haller
 * \version 0.1
 * \date    22.11.2007
 **************************************************************************/
#include "tom/core/volume.hpp"


#include <string.h>
#include <limits>
#include <iostream>
#include <assert.h>
#include <fftw3.h>


#include "tom/core/volume_loop.hpp"
#include "tom/core/transform.hpp"
#include "tom/core/volume_fcn.hpp"









namespace tom {
/****************************************************************************//**
 * \brief Structure where every tom::Volume remembers whether to free its memory.
 *
 * Volume can share the same memory using the constructor. All these objects point
 * to the same st_volume_memory to know how many other objects shares the memory
 * with them, and whether and how to free it.
 *******************************************************************************/
struct st_volume_memory {
    void *ptr;                      /**< Pointer to the memory which can be passed to free. Must not be the same as the base pointer of its volume. */
    std::size_t size;               /**< The maximum size allocated in bytes, beginning from ptr. */
    bool free;                      /**< If true, the destructor of the last object frees the memory calling fcn_free. */
    int cnt;                        /**< The number of objects sharing the same memory. */
    void (*fcn_free)(void *ptr);    /**< A pointer to the deallocation function. if !free it is NULL, otherwise it must be set. */
};
}


/****************************************************************************//**
 * \brief Returns a structure containing the memory infos of the volume.
 *******************************************************************************/
template<typename T>
std::size_t tom::Volume<T>::get_memsize() const {
    return this->volume_memory->size - ((const char *)this->data - (const char *)this->volume_memory->ptr);
}





/* saves default malloc-function.
 * Returned by tom::getFcnMalloc() */
void *(*fcn_malloc_default)(size_t n ) = NULL;

/* saves default free-function.
 * Returned by tom::getFcnMalloc() */
void  (*fcn_free_default  )(void *ptr) = NULL;






/****************************************************************************//**
 * \brief Sets the default memory allocation functions.
 *
 * \param[in] fcn_malloc Function to allocate memory. Synopsis must be the same
 *   as for malloc() from ISO C89.
 * \param[in] fcn_free Function to de-allocate memory. Synopsis must be the same
 *   as for free() from ISO C89.
 *
 * Sets the allocation and de-allocation functions as used in some cases to
 * allocate memory on the heap. The user must take care, that they correspond,
 * that means, memory allocated by \c fcn_malloc can be freed calling \c fcn_free and
 * vice versa.
 * Both parameters can be set to \c NULL. In that case the methods from tom_volume.cpp
 * use the new[] and delete[] operaters from C++. Not that these call
 * the default constructor of T, while malloc only allocates the memory.\n
 * Upon program start these functions are initialized to \c NULL.
 * In case of failure malloc should return NULL or throw std::bad_alloc. free should
 * fail or throw exceptions.
 *******************************************************************************/
void tom::setFcnMem( void *(*fcn_malloc)(size_t n),
                void (*fcn_free)(void *ptr)) {
    if ((!fcn_malloc && fcn_free) || (fcn_malloc && !fcn_free)) {
        throw std::invalid_argument("You must specify both fcn_malloc() and fcn_free() or leave both unspecified.");
    }
    fcn_malloc_default = fcn_malloc;
    fcn_free_default = fcn_free;
}



/****************************************************************************//**
 * \brief Returns the default memory allocation function.
 *
 * Can be set by tom::setMemFcn() (even to NULL).
 *******************************************************************************/
void *(*(tom::getFcnMalloc()))(size_t n) {
    return fcn_malloc_default;
}

/****************************************************************************//**
 * \brief Returns the default memory de-allocation function.
 *
 * Can be set by tom::setMemFcn() (even to NULL).
 *******************************************************************************/
void (*(tom::getFcnFree()))(void *) {
    return fcn_free_default;
}




/****************************************************************************//**
 * \brief Returns true, if the objects destructor will free the memory.
 *
 * Of course the memory is only freed it there are not other volumes left which
 * use the same data. See also number_memory_used().
 *******************************************************************************/
template<typename T>
bool tom::Volume<T>::to_free_memory() const {
    return this->volume_memory->free;
}


/****************************************************************************//**
 * \brief Returns how many objects in total share this data.
 *
 * This value will always be larger or equal to 1. The destructor of the object
 * decrements this value. If to_free_memory() is true and number_memory_used()
 * is 0, the destructor calls the de-allocation function.
 *******************************************************************************/
template<typename T>
int tom::Volume<T>::number_memory_used() const {
    return this->volume_memory->cnt;
}


/****************************************************************************//**
 * \brief Destructor of the volume.
 *
 * The destructor of the object the number of objects sharing the same data.
 * If to_free_memory() is true and number_memory_used()
 * is 0, the destructor calls the de-allocation function.
 *******************************************************************************/
template<typename T>
tom::Volume<T>::~Volume() {
    //std::cout << "call destructor: "; this->printInfo("v(destr)");
    if (--this->volume_memory->cnt < 1) {
        if (this->volume_memory->free) {
            //std::cout << "  +++++ do free memory" << std::endl;
            this->volume_memory->fcn_free(this->volume_memory->ptr);
        }
        //std::cout << "  +++++ free volume_memory" << std::endl;
        delete this->volume_memory;
    }
}







/****************************************************************************//**
 * \brief Initializes the volume by taking the data from an existing one.
 *
 * \param[in] v A reference to the volume which memory will be taken.
 * \param[in] data A pointer to the first element of the new volume. Setting
 *   this parameter to NULL, means taking the first element of v. (= &v.get()).
 * \param[in] sizex Size of the new volume.
 * \param[in] sizey Size of the new volume.
 * \param[in] sizez Size of the new volume.
 * \param[in] stridex Stride along x (in bytes).
 * \param[in] stridey Stride along y (in bytes).
 * \param[in] stridez Stride along z (in bytes).
 *******************************************************************************/
template<typename T> template<typename T2>
void tom::Volume<T>::share_memory(Volume<T2> &v, void *data, std::size_t sizex, std::size_t sizey, std::size_t sizez, std::size_t stridex, std::size_t stridey, std::size_t stridez) {

    if (sizex<=0 || sizey<=0 || sizez<=0) {
        throw std::invalid_argument("tom::Volume::share_memory - Empty volume is not allowed.");
    }
    if (!tom_io_check_stride(sizeof(T), sizex, sizey, sizez, &stridex, &stridey, &stridez)) {
        throw std::invalid_argument("tom::Volume::share_memory - Stride paramter to small according to data-type and volume size.");
    }

    std::size_t i;

    const char *database = reinterpret_cast<char *>(v.volume_memory->ptr);
    char *datac = reinterpret_cast<char *>(data);
    if (!datac) {
        datac = reinterpret_cast<char *>(&v.get());
    }

    if ((  datac < database) ||
        (i=datac - database) > v.volume_memory->size ||
        (v.volume_memory->size - i) < ((sizez-1) *stridez + (sizey-1)*stridey + (sizex-1) * stridex + sizeof(T))) {
        throw std::out_of_range("tom::Volume::share_memory - New data lies outside the (known) allocated memory.");
    }


    this->volume_memory = v.volume_memory;
    this->volume_memory->cnt++;

    this->data = const_cast<T *>(reinterpret_cast<const T *>(datac));
    this->sizex = sizex;
    this->sizey = sizey;
    this->sizez = sizez;
    this->stridex = stridex;
    this->stridey = stridey;
    this->stridez = stridez;

}


/****************************************************************************//**
 * \brief Create an empty volume.
 * This constructor is private, since no empty volumes are allowed.
 * This constructor can be used, to create an non initialized volume within
 * a member function or friend.
 *******************************************************************************/
template<typename T>
tom::Volume<T>::Volume()
    :   data(NULL),
        sizex(0),
        sizey(0),
        sizez(0),
        stridex(0),
        stridey(0),
        stridez(0),
        volume_memory(NULL) {
}


/****************************************************************************//**
 * \brief Creates a volume from already allocated memory.
 *
 * \param[in] data A pointer to the volume.
 * \param[in] sizex The size of the volume. Must be positive.
 * \param[in] sizey The size of the volume. Must be positive.
 * \param[in] sizez The size of the volume. Must be positive.
 * \param[in] free Whether the last volume should free the allocated memory in its
 *   destructor.
 * \param[in] stridex The distance in memory between two x-elements. This is
 *   measured in bytes, thus for contiguous memory, this equals to sizeof(T).
 *   Setting to zero defaults to contiguous
 * \param[in] stridey The distance in memory between two neighoubing y-elements.
 *   in bytes. For contiguous memory, this equals to sizex*stridex.
 * \param[in] stridez Same for z.
 * \param[in] fcn_free The function to de-allocate the memory afterwards.
 *   \c NULL defaults to tom::getFcnFree(), which by itself defaults to the
 *   C++ delete[] T operator, if set to \c NULL. In that case the function
 *   which is returned from tom::getFcnFree() at construction time is used,
 *   not the one from deconstruction time.
 *
 * You are responsable, that the memory pointed by data is large enough to
 * hold all values. Further the pointer must be valid as long as the volume
 * exists, and must be freed by the volume itself (in case of \c free ) or
 * by the caller.\n
 * \c fcn_free must be the right function to deallocate the memory. If the
 * memory was allocated using the new[] operator, set \c fcn_free to \c NULL
 * and tom::setFcnMem() also.\n
 * If you want to share memory with an other volume, use the appropriate
 * constructor instead.
 *******************************************************************************/
template<typename T>
tom::Volume<T>::Volume(T *data, std::size_t sizex, std::size_t sizey, std::size_t sizez, std::size_t stridex, std::size_t stridey, std::size_t stridez, bool free, void (*fcn_free)(void *ptr)) {

    if (!data || sizex<=0 || sizey<=0 || sizez<=0) {
        throw std::invalid_argument("Empty volume is not allowed.");
    }

    if (!tom_io_check_stride(sizeof(T), sizex, sizey, sizez, &stridex, &stridey, &stridez)) {
        throw std::invalid_argument("Stride paramter to small according to data-type and volume size.");
    }

    if (free) {
        if (!fcn_free) {
            if (!(fcn_free = tom::getFcnFree())) {
                fcn_free = &tom::fcn_free_delete<T>;
            }
        }
    } else {
        fcn_free = NULL;
    }

    this->volume_memory = new st_volume_memory();
    this->volume_memory->ptr = data;
    this->volume_memory->size = (sizez-1) *stridez + (sizey-1)*stridey + sizex * stridex;
    this->volume_memory->free = free;
    this->volume_memory->cnt = 1;
    this->volume_memory->fcn_free = fcn_free;

    this->data = data;
    this->sizex = sizex;
    this->sizey = sizey;
    this->sizez = sizez;
    this->stridex = stridex;
    this->stridey = stridey;
    this->stridez = stridez;
}




/****************************************************************************//**
 * \brief Initialize volume by allocating new (continous) memory.
 *******************************************************************************/
template<typename T>
void tom::Volume<T>::initialize(std::size_t sizex, std::size_t sizey, std::size_t sizez, void *(*fcn_malloc)(size_t size), void (*fcn_free)(void *ptr)) {

    if (sizex<=0 || sizey<=0 || sizez<=0) {
        throw std::invalid_argument("Empty volume is not allowed.");
    }

    if (!fcn_malloc && !fcn_free) {
        fcn_malloc = tom::getFcnMalloc();
        fcn_free = tom::getFcnFree();
        if (!fcn_malloc || !fcn_free) {
            fcn_malloc = &tom::fcn_malloc_new<T>;
            fcn_free = &tom::fcn_free_delete<T>;
        }
    } else if (!(fcn_malloc && fcn_free)) {
        throw std::invalid_argument("You must specify both fcn_malloc() and fcn_free() or leave both unspecified.");
    }

    size_t numel_bytes = sizex*sizey*sizez*sizeof(T);

    try {
        data = (T *)fcn_malloc(numel_bytes);
    } catch (...) {
        throw std::bad_alloc();
    }
    if (!data) {
        throw std::bad_alloc();
    }

    try {
        this->volume_memory = new st_volume_memory();
    } catch (...) {
        fcn_free(data);
        throw; /* Rethrow exception. */
    }

    this->volume_memory->ptr = data;
    this->volume_memory->size = numel_bytes;
    this->volume_memory->free = true;
    this->volume_memory->cnt = 1;
    this->volume_memory->fcn_free = fcn_free;

    this->data = data;
    this->sizex = sizex;
    this->sizey = sizey;
    this->sizez = sizez;
    this->stridex = sizeof(T);
    this->stridey = this->sizex * this->stridex;
    this->stridez = this->sizey * this->stridey;

}


/****************************************************************************//**
 * \brief Writes the volume as em-file.
 *
 * \param[in] filename The name of the em-file.
 * \param[in] header A pointer to the used em-header. If set to NULL, a default
 *   header is used with only the needed fields set according to the calling
 *   object. Even if a header is given, the fields \c dims and \c type are set
 *   to match the data from the volume.
 *
 * Calles tom_io_em_write to write the data. If an error occures
 * there, its integer error status is thrown as exception.
 *******************************************************************************/
template<typename T>
void tom::Volume<T>::write_to_em(const std::string &filename, const tom_io_em_header *header) const {

    tom_io_em_header header_local;

    if (this->getSizeX() > std::numeric_limits<uint32_t>::max() ||
        this->getSizeY() > std::numeric_limits<uint32_t>::max() ||
        this->getSizeZ() > std::numeric_limits<uint32_t>::max()) {
        throw std::runtime_error("The volume is too large to save it as em-format.");
    }


    if (header) {
        header_local = *header;
    } else {
        memset(&header_local, 0, sizeof(header_local));
        header_local.machine = 6;
    }

    if (!tom_io_em_set_iotype(&header_local, tom::get_tom_io_type<T>())) {
        throw std::runtime_error("This type can not be saved to em-file.");
    }
    header_local.dims[0] = this->getSizeX();
    header_local.dims[1] = this->getSizeY();
    header_local.dims[2] = this->getSizeZ();

    int res  = tom_io_em_write(filename.c_str(), &header_local, this->data, tom_io_em_get_iotype(&header_local), this->getStrideX(), this->getStrideY(), this->getStrideZ());

    if (res != TOM_ERR_OK) {
        throw res;
    }
}

/****************************************************************************//**
 * \brief Writes the volume as em-file.
 *
 * \param[in] filename The name of the em-file.
 * \param[in] header A pointer to the used em-header. If set to NULL, a default
 *   header is used with only the needed fields set according to the calling
 *   object. Even if a header is given, the fields \c dims and \c type are set
 *   to match the data from the volume.
 * \param[in] first_voxel position where the volume is written to.
 *
 * Calles tom_io_em_write to write the data. If an error occures
 * there, its integer error status is thrown as exception.
 *******************************************************************************/
template<typename T>
void tom::Volume<T>::write_to_em(const std::string &filename, const tom_io_em_header *header,const uint32_t* first_voxel) const {

    tom_io_em_header header_local;

    if (this->getSizeX() > std::numeric_limits<uint32_t>::max() ||
        this->getSizeY() > std::numeric_limits<uint32_t>::max() ||
        this->getSizeZ() > std::numeric_limits<uint32_t>::max()) {
        throw std::runtime_error("tom::Volume<T>::write_to_em-paste : The volume is too large to save it as em-format.");
    }


    if (header) {
        header_local = *header;
    } else {
        memset(&header_local, 0, sizeof(header_local));
        header_local.machine = 6;
    }

    if (!tom_io_em_set_iotype(&header_local, tom::get_tom_io_type<T>())) {
        throw std::runtime_error("tom::Volume<T>::write_to_em-paste : This type can not be saved to em-file.");
    }
    header_local.dims[0] = this->getSizeX();
    header_local.dims[1] = this->getSizeY();
    header_local.dims[2] = this->getSizeZ();

	int header_read;

	int res  = tom_io_em_write_paste(filename.c_str(),this->data,tom_io_em_get_iotype(&header_local),
									this->getSizeX(), this->getSizeY(), this->getSizeZ(),
									this->getStrideX(), this->getStrideY(),this->getStrideZ(),
									&header_local, &header_read, true, first_voxel, NULL);

    if (res != TOM_ERR_OK) {
        throw res;
    }
}


/****************************************************************************//**
 * \brief A malloc function using internally the C++ new[] operator.
 *
 * \param[in] size The number of bytes to be allocated.
 *
 * Creates <tt> ceil(size/sizeof(T)) </tt> objects of \c T with its default
 * constructor. Which should make not a big difference for elementary datatypes.
 *******************************************************************************/
template<typename T> void *tom::fcn_malloc_new(size_t size) {
    size_t remainder = size % sizeof(T);
    size /= sizeof(T);
    if (remainder > 0) {
        size++;
    }
    return new T[size];
}

/****************************************************************************//**
 * \brief A free function using internally the C++ delete[] operator.
 *******************************************************************************/
template<typename T> void tom::fcn_free_delete(void *ptr) {
    delete[] (T *)ptr;
}


/****************************************************************************//**
 * \brief Reads a volume from em-file.
 *
 * \param[out] v A reference to a pointer where the volume is returned.
 * \param[in] filename The filename of the emfile.
 * \param[in] subregion The subregion to read. See tom_io_em_read.
 * \param[in] sampling The sampling factor. See tom_io_em_read.
 * \param[in] binning The binning factor. See tom_io_em_read.
 * \param[out] header A pointer to the em-header. Set to \c NULL if you are not
 *   interested in the header.
 * \param[in] fcn_malloc A pointer to a function with which the memory will be
 *   allocated. Setting to \c NULL defaults to tom::getFcnMalloc(). If
 *   getFcnMalloc returns \c NULL too, the new[] for type T is used.
 * \param[in] fcn_free A pointer to a function with which the memory allocated with
 *   fcn_malloc can be freed again. Setting to \c NULL defaults to
 *   tom::getFcnFree(), and here too, to delete[].
 *
 * Calls tom_io_em_read to read the file, and throws the integer error status
 * returned by it, in case of an error.
 *******************************************************************************/
template<typename T>
void tom::read_from_em(Volume<T> *&v, const std::string &filename, const uint32_t *subregion, const uint32_t *sampling, const uint32_t *binning, tom_io_em_header *header, void *(*fcn_malloc)(size_t size), void (*fcn_free)(void *ptr)) {

    if (!fcn_malloc && !fcn_free) {
        fcn_malloc = tom::getFcnMalloc();
        fcn_free = tom::getFcnFree();
        if (!fcn_malloc || !fcn_free) {
            fcn_malloc = &tom::fcn_malloc_new<T>;
            fcn_free = &tom::fcn_free_delete<T>;
        }
    } else if (!(fcn_malloc && fcn_free)) {
        throw std::invalid_argument("You must specify both fcn_malloc() and fcn_free() or leave both unspecified.");
    }


    void *vdata = NULL;
    uint32_t dims[3];
    int restype;
    std::auto_ptr<Volume<T> > v_local;

    try {
        int res = ::tom_io_em_read(filename.c_str(), subregion, sampling, binning, header, &vdata, dims, &restype, fcn_malloc, fcn_free);

        if (res != TOM_ERR_OK) {
            throw res;
        }

        if (restype == tom::get_tom_io_type<T>()) {
            v_local.reset(new tom::Volume<T>((T *)vdata, dims[0], dims[1], dims[2], 0, 0, 0, true, fcn_free));
            vdata = NULL;
        } else {
            v_local.reset(new tom::Volume<T>(dims[0], dims[1], dims[2], fcn_malloc, fcn_free));
            switch (restype) {
                case TOM_IO_TYPE_FLOAT: {
                        Volume<float> v2((float *)vdata, dims[0], dims[1], dims[2], 0, 0, 0, true, fcn_free);
                        vdata = NULL; /* The data from vdata is now owned by v2. Leave clean up to local object. */
                        if (tom::is_double<T>()) {
                            tom::Volume<double> *v_local_typed = reinterpret_cast<tom::Volume<double> *>(v_local.get());
                            v_local_typed->setValues(v2);
                        } else {
                            throw std::runtime_error("This type conversion is not (yet?) not implemented.");
                        }
                    }
                    break;
                case TOM_IO_TYPE_DOUBLE: {
                        Volume<double> v2((double *)vdata, dims[0], dims[1], dims[2], 0, 0, 0, true, fcn_free);
                        vdata = NULL; /* The data from vdata is now owned by v2. Leave it to him to clean up. */
                        if (tom::is_float<T>()) {
                            tom::Volume<float> *v_local_typed = reinterpret_cast<tom::Volume<float> *>(v_local.get());
                            v_local_typed->setValues(v2);
                        } else {
                            throw std::runtime_error("This type conversion is not (yet?) not implemented.");
                        }
                    }
                    break;
                case TOM_IO_TYPE_INT8:
                case TOM_IO_TYPE_INT16:
                case TOM_IO_TYPE_INT32:
                case TOM_IO_TYPE_COMPLEX4:
                default:
                    throw std::logic_error("This type conversion is not (yet?) not implemented.");
            }
        }
    } catch (int &a) {
        /* Clean up... */
        std::stringstream ss;
        ss << "Error reading the em-file \"" << filename << "\" (errorcode " << a << ").";
        if  (vdata) { fcn_free(vdata); }
        throw std::runtime_error(ss.str());
    } catch (...) {
        /* Clean up... */
        if  (vdata) { fcn_free(vdata); }
        throw;
    }
    v = v_local.release();
}






/****************************************************************************//**
 * \brief Writes debugging info to stdout.
 *******************************************************************************/
template<typename T>
void tom::Volume<T>::printInfo(const std::string &name) const {
    std::cout << "Volume " << name << ":" << std::endl <<
                 "  size[zyx] = [" << this->getSizeZ()   << "," << this->getSizeY()   << "," << this->getSizeX()   << "]; of tom_io_type " << tom::get_tom_io_type<T>() << std::endl <<
                 "  stride[zyx] = [" << this->getStrideZ() << "," << this->getStrideY() << "," << this->getStrideX() << "]; gap[zyx] = [" << (this->getStrideZ()-this->getStrideY()*this->getSizeY()) << "," << (this->getStrideY()-this->getStrideX()*this->getSizeX()) << "," << (this->getStrideX()-sizeof(T)) << "]; " << std::endl;
    if (tom::is_double<T>() || tom::is_float<T>()) {
        double mean, variance;
        T min, max;
        this->stat(mean, variance, min, max, true);
        std::cout << "  range = [" << min << " .. " << mean << " .. " << max << "]" << std::endl
                  << "  stddev = " << sqrt(this->variance(false)) << " (stddev_sample = " << sqrt(variance) << ")" << std::endl;
    }
    std::cout << "  ptr = " << &this->get() << ", baseptr = " << this->volume_memory->ptr << " (offset=" << (((size_t)&this->get()) - ((size_t)(this->volume_memory->ptr))) << ")" << std::endl <<
                 "  pointed to and ";
    if (!this->to_free_memory()) { std::cout << "not "; }
    std::cout << "to free by " << this->number_memory_used() << " volumes using " << ((void *)((size_t)this->volume_memory->fcn_free)) << " (";
    void (*fcn_free_delete_double)(void*) = &tom::fcn_free_delete<double>;
    void (*fcn_free_delete_float)(void*) = &tom::fcn_free_delete<float>;
    if (this->volume_memory->fcn_free == &fftw_free) {
        std::cout << "fftw_free";
    } else if (this->volume_memory->fcn_free == &fftwf_free) {
        std::cout << "fftwf_free";
    } else if (this->volume_memory->fcn_free == fcn_free_delete_double) {
        std::cout << "fcn_free_delete<double>";
    } else if (this->volume_memory->fcn_free == fcn_free_delete_float) {
        std::cout << "fcn_free_delete<float>";
    } else if (this->volume_memory->fcn_free == &free) {
        std::cout << "free";
    } else if (!this->volume_memory->fcn_free) {
        std::cout << "NULL";
    } else {
        std::cout << "unknown";
    }
    std::cout << ")" << std::endl;
}






namespace {
template<typename T> inline void setValues__contiguous_copy(T                    *dest, const T                    *src, std::size_t n) { for (std::size_t i=0; i<n; i++) { dest[i] = src[i]; } } // Dont use memcpy for calling the assignment operator.
template<          > inline void setValues__contiguous_copy(int                  *dest, const int                  *src, std::size_t n) { memcpy(dest, src, n*sizeof(int                 )); }
template<          > inline void setValues__contiguous_copy(float                *dest, const float                *src, std::size_t n) { memcpy(dest, src, n*sizeof(float               )); }
template<          > inline void setValues__contiguous_copy(double               *dest, const double               *src, std::size_t n) { memcpy(dest, src, n*sizeof(double              )); }
template<          > inline void setValues__contiguous_copy(std::complex<float > *dest, const std::complex<float > *src, std::size_t n) { memcpy(dest, src, n*sizeof(std::complex<float >)); }
template<          > inline void setValues__contiguous_copy(std::complex<double> *dest, const std::complex<double> *src, std::size_t n) { memcpy(dest, src, n*sizeof(std::complex<double>)); }
template<typename T1, typename T2> inline void for_each__tom__volume__setValues2_(T1 &v1, const T2 &v2) { v1 = v2; }
template<typename T1, typename T2>
struct for_each__tom__volume__setValues2 {
    inline void operator()(T1 &v1, const T2 &v2, std::size_t, std::size_t, std::size_t) { for_each__tom__volume__setValues2_<T1,T2>(v1, v2); }
};
}
/****************************************************************************//**
 * \brief Copy the values from one Volume to an other with type conversion.
 *
 * \param[in] v The source volume, of the same size as *this. The datatypes can
 *   differ.
 *******************************************************************************/
template<typename T>
template<typename T2>
void tom::Volume<T>::setValues(const Volume<T2> &v) {
    if (!this->is_equal_size(v)) {
        throw std::invalid_argument("Both volumes must have the same size.");
    }
    if (static_cast<const void *>(&v.get()) == static_cast<void *>(&this->get())) {
        if (typeid(T) != typeid(T2)) {
            throw std::invalid_argument("Self assigment of different datatypes.");
        }
        return;
    }
    tom::loop::for_each<T, tom::Volume<T> &, const T2, const tom::Volume<T2> &, ::for_each__tom__volume__setValues2<T, T2> >(*this, v, ::for_each__tom__volume__setValues2<T, T2>());
}
/****************************************************************************//**
 * \brief Copy the values from one Volume to an other of the same type.
 *
 * \param[in] v The source volume, of the same size as *this.
 *******************************************************************************/
template<typename T>
void tom::Volume<T>::setValues(const Volume<T> &v) {
    if (!this->is_equal_size<T>(v)) {
        throw std::invalid_argument("Both volumes must have the same size.");
    }
    if (&v.get() == &this->get()) {
        return;
    }

    if (this->isContiguous() && v.isContiguous()) {
        ::setValues__contiguous_copy<T>(&this->get(), &v.get(), this->getSizeX()*this->getSizeY()*this->getSizeZ());
    } else {
        tom::loop::for_each<T, tom::Volume<T> &, const T, const tom::Volume<T> &, ::for_each__tom__volume__setValues2<T, T> >(*this, v, ::for_each__tom__volume__setValues2<T, T>());
    }
}






namespace {
template<typename T>
struct for_each__tom__volume__setValues {
    for_each__tom__volume__setValues(T val): val(val) { }
    T val;
    inline void operator()(T &a, std::size_t, std::size_t, std::size_t) { a = this->val; }
};
}
/****************************************************************************//**
 * \brief Set all voxels to a specific value.
 *
 * \param[in] val The value which is assigned to each voxel of this.
 *******************************************************************************/
template<typename T>
void tom::Volume<T>::setValues(T val) {
    tom::loop::for_each<T, tom::Volume<T> &, ::for_each__tom__volume__setValues<T> >(*this, ::for_each__tom__volume__setValues<T>(val));
}




namespace {
template<typename T>
struct for_each_while__tom__volume__opequal {
    for_each_while__tom__volume__opequal(T val): val(val), equal(true) { }
    T val;
    bool equal;
    inline bool operator()(const T &a, std::size_t, std::size_t, std::size_t) { if (a != this->val) { this->equal = false; } return this->equal; }
};
}
/****************************************************************************//**
 * \brief Checks wether all elements of the volume equal one value.
 *
 * \param[in] val The value to be compared with the volume.
 * \returns true if all elements of the volume are equal the scalar.
 *   Otherwise false.
 *
 * If the first element not equal val is found, false is returned. otherwise
 * every value has to be checked.
 *******************************************************************************/
template<typename T>
bool tom::Volume<T>::operator==(const T &val) const {
    ::for_each_while__tom__volume__opequal<T> s(val);
    tom::loop::for_each_while<const T, const tom::Volume<T> &, ::for_each_while__tom__volume__opequal<T> &>(*this, s);
    return s.equal;
}




namespace {
template<typename T>
struct for_each_while__tom__volume__opequalv {
    for_each_while__tom__volume__opequalv(): equal(true) { }
    bool equal;
    inline bool operator()(const T &v1, const T &v2, std::size_t, std::size_t, std::size_t) { if (v1 != v2) { this->equal = false; } return this->equal; }
};
}
/****************************************************************************//**
 * \brief Checks wether all elements are equal to each other.
 *
 * \param[in] vol Volume to be compared.
 * \returns true if all elements of the volume equal.
 *   Otherwise false. If the volumes have different size, false is returned.
 *
 * If the first element not equal val is found, false is returned. otherwise
 * every value has to be checked.
 *******************************************************************************/
template<typename T>
bool tom::Volume<T>::operator==(const tom::Volume<T> &v) const {
    if (!this->is_equal_size<T>(v)) {
        return false;
    }
    if (this==&v || (&this->get()==&v.get() && this->getStrideX()==v.getStrideX() && this->getStrideY()==v.getStrideY() && this->getStrideZ()==v.getStrideZ())) {
        return true;
    }
    ::for_each_while__tom__volume__opequalv<T> s;
    tom::loop::for_each_while<const T, const tom::Volume<T> &, const T, const tom::Volume<T> &, ::for_each_while__tom__volume__opequalv<T> &>(*this, v, s);

    return s.equal;
}






namespace {
template<typename T, typename TPRECISION>
struct for_each__tom__volume__stat__minmaxsum2 {
    for_each__tom__volume__stat__minmaxsum2(T init): min(init), max(init), sum(0), sum2(0) { }
    T min, max;
    TPRECISION sum, sum2;
    inline void operator()(const T &a, std::size_t, std::size_t, std::size_t) {
        if (a < this->min) {
            this->min = a;
        } else if (a > this->max) {
            this->max = a;
        }
        this->sum += a;
        this->sum2 += a*a;
    }
};
template<typename T, typename TPRECISION>
struct for_each__tom__volume__stat__minmaxsum {
    for_each__tom__volume__stat__minmaxsum(T init): min(init), max(init), sum(0) { }
    T min, max;
    TPRECISION sum;
    inline void operator()(const T &a, std::size_t, std::size_t, std::size_t) {
        if (a < this->min) {
            this->min = a;
        } else if (a > this->max) {
            this->max = a;
        }
        this->sum += a;
    }
};
template<typename T>
struct for_each__tom__volume__stat__minmax {
    for_each__tom__volume__stat__minmax(T init): min(init), max(init) { }
    T min, max;
    inline void operator()(const T &a, std::size_t, std::size_t, std::size_t) {
        if (a < this->min) {
            this->min = a;
        } else if (a > this->max) {
            this->max = a;
        }
    }
};
template<typename T>
struct for_each__tom__volume__stat__min {
    for_each__tom__volume__stat__min(T init): min(init) { }
    T min;
    inline void operator()(const T &a, std::size_t, std::size_t, std::size_t) { if (a < this->min) { this->min = a; } }
};
template<typename T>
struct for_each__tom__volume__stat__max {
    for_each__tom__volume__stat__max(T init): max(init) { }
    T max;
    inline void operator()(const T &a, std::size_t, std::size_t, std::size_t) { if (a > this->max) { this->max = a; } }
};
template<typename T, typename TPRECISION>
struct for_each__tom__volume__stat__sum {
    for_each__tom__volume__stat__sum(): sum(0) { }
    TPRECISION sum;
    inline void operator()(const T &a, std::size_t, std::size_t, std::size_t) { this->sum += a; }
};
template<typename T, typename TPRECISION>
struct for_each__tom__volume__stat__sum2 {
    for_each__tom__volume__stat__sum2(): sum(0), sum2(0) { }
    TPRECISION sum, sum2;
    inline void operator()(const T &a, std::size_t, std::size_t, std::size_t) { this->sum += a; this->sum2 += a*a; }
};
}
/****************************************************************************//**
 * \brief Private method which computes the mean, variance and min/max.
 *
 * Depending on which parameter are needed, different things are calculated.
 *******************************************************************************/
template<typename T>
void tom::Volume<T>::stat(double *m, double *variance, T *min, T *max, bool use_sample_standard_deviation) const {

    typedef double TPRECISION;

    const std::size_t numel = this->getSizeX() * this->getSizeY() * this->getSizeZ();

    TPRECISION lmean = 0.;
    TPRECISION lvariance = 0.;
    T lmin, lmax;
    lmin = lmax = this->get();

    if (m && !variance && !min && !max) {
        ::for_each__tom__volume__stat__sum<T, TPRECISION> s;
        tom::loop::for_each<const T, const tom::Volume<T> &, ::for_each__tom__volume__stat__sum<T, TPRECISION> &>(*this, s);
        lmean = s.sum / static_cast<TPRECISION>(numel);
    } else if (variance && !min && !max) {
        ::for_each__tom__volume__stat__sum2<T, TPRECISION> s;
        tom::loop::for_each<const T, const tom::Volume<T> &, ::for_each__tom__volume__stat__sum2<T, TPRECISION> &>(*this, s);
        lmean = s.sum / static_cast<TPRECISION>(numel);
        lvariance = (s.sum2 - static_cast<TPRECISION>(numel)*lmean*lmean) / (static_cast<TPRECISION>(numel)- (use_sample_standard_deviation ? 1 : 0));
    } else if (!m && !variance && min && !max) {
        ::for_each__tom__volume__stat__min<T> s(this->get());
        tom::loop::for_each<const T, const tom::Volume<T> &, ::for_each__tom__volume__stat__min<T> &>(*this, s);
        lmin = s.min;
    } else if (!m && !variance && !min && max) {
        ::for_each__tom__volume__stat__max<T> s(this->get());
        tom::loop::for_each<const T, const tom::Volume<T> &, ::for_each__tom__volume__stat__max<T> &>(*this, s);
        lmax = s.max;
    } else if (!m && !variance && min && max) {
        ::for_each__tom__volume__stat__minmax<T> s(this->get());
        tom::loop::for_each<const T, const tom::Volume<T> &, ::for_each__tom__volume__stat__minmax<T> &>(*this, s);
        lmin = s.min;
        lmax = s.max;
    } else if (!variance) {
        ::for_each__tom__volume__stat__minmaxsum<T, TPRECISION> s(this->get());
        tom::loop::for_each<const T, const tom::Volume<T> &, ::for_each__tom__volume__stat__minmaxsum<T, TPRECISION> &>(*this, s);
        lmin = s.min;
        lmax = s.max;
        lmean = s.sum / static_cast<TPRECISION>(numel);
    } else if (m || variance || min || max) {
        ::for_each__tom__volume__stat__minmaxsum2<T, TPRECISION> s(this->get());
        tom::loop::for_each<const T, const tom::Volume<T> &, ::for_each__tom__volume__stat__minmaxsum2<T, TPRECISION> &>(*this, s);
        lmin = s.min;
        lmax = s.max;
        lmean = s.sum / static_cast<TPRECISION>(numel);
        lvariance = (s.sum2 - static_cast<TPRECISION>(numel)*lmean*lmean) / (static_cast<TPRECISION>(numel)- (use_sample_standard_deviation ? 1 : 0));
    }
    if (max) { *max = lmax; }
    if (min) { *min = lmin; }
    if (m) { *m = lmean; }
    if (variance) { *variance = lvariance; }
}
/** \cond __HIDE_FUNCTIONS */
namespace tom {
template<> void Volume<std::complex<float > >::stat(double *m, double *variance, std::complex<float > *min, std::complex<float > *max, bool use_sample_standard_deviation) const { throw std::runtime_error("function stat not defined for complex volume"); }
template<> void Volume<std::complex<double> >::stat(double *m, double *variance, std::complex<double> *min, std::complex<double> *max, bool use_sample_standard_deviation) const { throw std::runtime_error("function stat not defined for complex volume"); }
}
/** \endcond __HIDE_FUNCTIONS */



namespace {
template<typename T, typename TPRECISION>
struct for_each__tom__volume__variance_mean_free {
    for_each__tom__volume__variance_mean_free(): sum2(0) { }
    TPRECISION sum2;
    inline void operator()(const T &a, std::size_t, std::size_t, std::size_t) { this->sum2 += a*a; }
};
}
/****************************************************************************//**
 * \brief Computes the variance of a volume which has a mean of zero.
 *
 * \param[in] use_sample_standard_deviation If true, the variance is computed
 *    using the standard deviation of the sample (with nominator N-1).
 *
 * \return The variance of the volume.
 *
 * If the mean of the volume is not 0, the result is wrong. This case is not
 * checked.
 *******************************************************************************/
template<typename T>
double tom::Volume<T>::variance_mean_free(bool use_sample_standard_deviation) const {

    typedef double TPRECISION;

    ::for_each__tom__volume__variance_mean_free<T, TPRECISION> s;
    tom::loop::for_each<const T, const tom::Volume<T> &, ::for_each__tom__volume__variance_mean_free<T, TPRECISION> &>(*this, s);

    return s.sum2 / static_cast<TPRECISION>(this->getSizeX() * this->getSizeY() * this->getSizeZ() - (use_sample_standard_deviation ? 1 : 0));
}
/** \cond __HIDE_FUNCTIONS */
namespace tom {
template<> double Volume<std::complex<float > >::variance_mean_free(bool use_sample_standard_deviation) const { throw std::runtime_error("function stat not defined for complex volume"); }
template<> double Volume<std::complex<double> >::variance_mean_free(bool use_sample_standard_deviation) const { throw std::runtime_error("function stat not defined for complex volume"); }
}
/** \endcond __HIDE_FUNCTIONS */



namespace {
template<typename T, typename TPRECISION>
inline void for_each__tom__volume__shift_scale_(T &a, const TPRECISION &factor_shift, const TPRECISION &factor_scale) {
    a = (a + factor_shift) * factor_scale;
}
template<> inline void for_each__tom__volume__shift_scale_<std::complex<float >, float >(std::complex<float > &a, const float  &factor_shift, const float  &factor_scale) { a = (a                       +                     factor_shift        ) *              factor_scale ; }
template<> inline void for_each__tom__volume__shift_scale_<std::complex<float >, double>(std::complex<float > &a, const double &factor_shift, const double &factor_scale) { a = (std::complex<double>(a) +                     factor_shift        ) *              factor_scale ; }
template<> inline void for_each__tom__volume__shift_scale_<std::complex<double>, float >(std::complex<double> &a, const float  &factor_shift, const float  &factor_scale) { a = (a                       + static_cast<double>(factor_shift)) * static_cast<double>(factor_scale); }
template<> inline void for_each__tom__volume__shift_scale_<std::complex<double>, double>(std::complex<double> &a, const double &factor_shift, const double &factor_scale) { a = (a                       +                     factor_shift        ) *              factor_scale ; }
template<typename T, typename TPRECISION>
struct for_each__tom__volume__shift_scale {
    for_each__tom__volume__shift_scale(TPRECISION factor_shift, TPRECISION factor_scale): factor_shift(factor_shift), factor_scale(factor_scale) { }
    TPRECISION factor_shift, factor_scale;
    inline void operator()(T &a, std::size_t, std::size_t, std::size_t) { ::for_each__tom__volume__shift_scale_<T, TPRECISION>(a, this->factor_shift, this->factor_scale); }
};
template<typename T, typename TPRECISION>
inline void for_each__tom__volume__shift_(T &a, const TPRECISION &factor_shift) {
    a = a + factor_shift;
}
template<> inline void for_each__tom__volume__shift_<std::complex<float >, float >(std::complex<float > &a, const float  &factor_shift) { a = a                       +                     factor_shift ; }
template<> inline void for_each__tom__volume__shift_<std::complex<float >, double>(std::complex<float > &a, const double &factor_shift) { a = std::complex<double>(a) +                     factor_shift ; }
template<> inline void for_each__tom__volume__shift_<std::complex<double>, float >(std::complex<double> &a, const float  &factor_shift) { a = a                       + static_cast<double>(factor_shift); }
template<> inline void for_each__tom__volume__shift_<std::complex<double>, double>(std::complex<double> &a, const double &factor_shift) { a = a                       +                     factor_shift ; }
template<typename T, typename TPRECISION>
struct for_each__tom__volume__shift {
    for_each__tom__volume__shift(TPRECISION factor): factor(factor) { }
    TPRECISION factor;
    inline void operator()(T &a, std::size_t, std::size_t, std::size_t) { ::for_each__tom__volume__shift_<T, TPRECISION>(a, this->factor); }
};

template<typename T, typename TPRECISION>
inline void for_each__tom__volume__scale_(T &a, const TPRECISION &factor_scale) {
    a = a * factor_scale;
}
template<> inline void for_each__tom__volume__scale_<std::complex<float >, float >(std::complex<float > &a, const float  &factor_scale) { a = a                       *                     factor_scale ; }
template<> inline void for_each__tom__volume__scale_<std::complex<float >, double>(std::complex<float > &a, const double &factor_scale) { a = std::complex<double>(a) *                     factor_scale ; }
template<> inline void for_each__tom__volume__scale_<std::complex<double>, float >(std::complex<double> &a, const float  &factor_scale) { a = a                       * static_cast<double>(factor_scale); }
template<> inline void for_each__tom__volume__scale_<std::complex<double>, double>(std::complex<double> &a, const double &factor_scale) { a = a                       *                     factor_scale ; }
template<typename T, typename TPRECISION>
struct for_each__tom__volume__scale {
    for_each__tom__volume__scale(TPRECISION factor): factor(factor) { }
    TPRECISION factor;
    inline void operator()(T &a, std::size_t, std::size_t, std::size_t) { ::for_each__tom__volume__scale_<T, TPRECISION>(a, this->factor); }
};
}
/****************************************************************************//**
 * \brief Adds and scales the volume by a scalar.
 *******************************************************************************/
template<typename T>
template<typename TPRECISION>
void tom::Volume<T>::shift_scale(TPRECISION factor_shift, TPRECISION factor_scale) {
    if (factor_shift!=0 && factor_scale!=1) {
        tom::loop::for_each<T, tom::Volume<T> &, ::for_each__tom__volume__shift_scale<T, TPRECISION> >(*this, ::for_each__tom__volume__shift_scale<T, TPRECISION>(factor_shift, factor_scale));
    } else if (factor_shift!=0) {
        tom::loop::for_each<T, tom::Volume<T> &, ::for_each__tom__volume__shift<T, TPRECISION> >(*this, ::for_each__tom__volume__shift<T, TPRECISION>(factor_shift));
    } else if (factor_scale!=1) {
        tom::loop::for_each<T, tom::Volume<T> &, ::for_each__tom__volume__scale<T, TPRECISION> >(*this, ::for_each__tom__volume__scale<T, TPRECISION>(factor_scale));
    }
}







// template instantiation.
/** \cond __HIDE_FUNCTIONS */

template class tom::Volume<char         >;
template class tom::Volume<int          >;
template class tom::Volume<float        >;
template class tom::Volume<double       >;
template class tom::Volume<std::complex<float > >;
template class tom::Volume<std::complex<double> >;

template tom::Volume<float                >::Volume(Volume<float                > &v, void *data, std::size_t sizex, std::size_t sizey, std::size_t sizez, std::size_t stridex, std::size_t stridey, std::size_t stridez);
template tom::Volume<float                >::Volume(Volume<std::complex<float > > &v, void *data, std::size_t sizex, std::size_t sizey, std::size_t sizez, std::size_t stridex, std::size_t stridey, std::size_t stridez);
template tom::Volume<double               >::Volume(Volume<double               > &v, void *data, std::size_t sizex, std::size_t sizey, std::size_t sizez, std::size_t stridex, std::size_t stridey, std::size_t stridez);
template tom::Volume<double               >::Volume(Volume<std::complex<double> > &v, void *data, std::size_t sizex, std::size_t sizey, std::size_t sizez, std::size_t stridex, std::size_t stridey, std::size_t stridez);
template tom::Volume<std::complex<float > >::Volume(Volume<std::complex<float > > &v, void *data, std::size_t sizex, std::size_t sizey, std::size_t sizez, std::size_t stridex, std::size_t stridey, std::size_t stridez);
template tom::Volume<std::complex<double> >::Volume(Volume<std::complex<double> > &v, void *data, std::size_t sizex, std::size_t sizey, std::size_t sizez, std::size_t stridex, std::size_t stridey, std::size_t stridez);



template void tom::Volume<float   >::setValues<char>(const Volume<char> &v);
template void tom::Volume<double  >::setValues<char>(const Volume<char> &v);
template void tom::Volume<float   >::setValues<int >(const Volume<int > &v);
template void tom::Volume<double  >::setValues<int >(const Volume<int > &v);



template void tom::Volume<float                >::shift_scale<float >(float  factor_shift, float  factor_scale);
template void tom::Volume<float                >::shift_scale<double>(double factor_shift, double factor_scale);
template void tom::Volume<double               >::shift_scale<double>(double factor_shift, double factor_scale);
template void tom::Volume<std::complex<float > >::shift_scale<float >(float  factor_shift, float  factor_scale);
template void tom::Volume<std::complex<float > >::shift_scale<double>(double factor_shift, double factor_scale);
template void tom::Volume<std::complex<double> >::shift_scale<float >(float  factor_shift, float  factor_scale);
template void tom::Volume<std::complex<double> >::shift_scale<double>(double factor_shift, double factor_scale);

template void tom::Volume<int   >::share_memory(Volume<int   > &v, void *data, std::size_t sizex, std::size_t sizey, std::size_t sizez, std::size_t stridex, std::size_t stridey, std::size_t stridez);
template void tom::Volume<double>::share_memory(Volume<float > &v, void *data, std::size_t sizex, std::size_t sizey, std::size_t sizez, std::size_t stridex, std::size_t stridey, std::size_t stridez);
template void tom::Volume<float >::share_memory(Volume<double> &v, void *data, std::size_t sizex, std::size_t sizey, std::size_t sizez, std::size_t stridex, std::size_t stridey, std::size_t stridez);


template void tom::read_from_em<int32_t>(Volume<int32_t> *&v, const std::string &filename, const uint32_t *subregion, const uint32_t *sampling, const uint32_t *binning, tom_io_em_header *header, void *(*fcn_malloc)(size_t size), void (*fcn_free)(void *ptr));
template void tom::read_from_em<float  >(Volume<float  > *&v, const std::string &filename, const uint32_t *subregion, const uint32_t *sampling, const uint32_t *binning, tom_io_em_header *header, void *(*fcn_malloc)(size_t size), void (*fcn_free)(void *ptr));
template void tom::read_from_em<double >(Volume<double > *&v, const std::string &filename, const uint32_t *subregion, const uint32_t *sampling, const uint32_t *binning, tom_io_em_header *header, void *(*fcn_malloc)(size_t size), void (*fcn_free)(void *ptr));

/** \endcond __HIDE_FUNCTIONS */

