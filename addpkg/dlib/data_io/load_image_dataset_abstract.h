// Copyright (C) 2012  Davis E. King (davis@dlib.net)
// License: Boost Software License   See LICENSE.txt for the full license.
#undef DLIB_LOAD_IMAGE_DaTASET_ABSTRACT_H__
#ifdef DLIB_LOAD_IMAGE_DaTASET_ABSTRACT_H__

#include "image_dataset_metadata.h"
#include "../array/array_kernel_abstract.h"
#include <string>
#include <vector>


namespace dlib
{

// ----------------------------------------------------------------------------------------

    template <
        typename image_type, 
        typename MM
        >
    void load_image_dataset (
        array<image_type,MM>& images,
        std::vector<std::vector<rectangle> >& object_locations,
        const std::string& filename,
        const std::string& label 
    );
    /*!
        requires
            - image_type == is an implementation of array2d/array2d_kernel_abstract.h
            - pixel_traits<typename image_type::type> is defined  
        ensures
            - This routine loads the images and their associated object boxes from 
              the image metadata file indicated by filename.  This metadata file
              should be in the XML format used by the save_image_dataset_metadata()
              routine.
            - #images.size() == the number of images in the metadata file
            - #images.size() == #object_locations.size()
            - This routine is capable of loading any image format which can be read
              by the load_image() routine.
            - for all valid i:  
                - #images[i] == a copy of the ith image from the dataset
                - #object_locations[i] == a vector of all the rectangles associated with
                  #images[i].
                - if (labels != "") then
                    - only boxes with the given label will be loaded into object_locations.
                - else
                    - all boxes in the dataset will be loaded into object_locations.
    !*/

// ----------------------------------------------------------------------------------------

    template <
        typename image_type, 
        typename MM
        >
    void load_image_dataset (
        array<image_type,MM>& images,
        std::vector<std::vector<rectangle> >& object_locations,
        const std::string& filename
    );
    /*!
        requires
            - image_type == is an implementation of array2d/array2d_kernel_abstract.h
            - pixel_traits<typename image_type::type> is defined  
        ensures
            - performs: load_image_dataset(images, object_locations, filename, "");
              (i.e. it ignores box labels and therefore loads all the boxes in the dataset)
    !*/

// ----------------------------------------------------------------------------------------

}

#endif // DLIB_LOAD_IMAGE_DaTASET_ABSTRACT_H__


