// Copyright (C) 2011  Davis E. King (davis@dlib.net)
// License: Boost Software License   See LICENSE.txt for the full license.
#undef DLIB_STRUCTURAL_OBJECT_DETECTION_TRAiNER_H_ABSTRACT__
#ifdef DLIB_STRUCTURAL_OBJECT_DETECTION_TRAiNER_H_ABSTRACT__

#include "structural_svm_object_detection_problem_abstract.h"
#include "../image_processing/object_detector_abstract.h"
#include "../image_processing/box_overlap_testing_abstract.h"


namespace dlib
{

// ----------------------------------------------------------------------------------------

    template <
        typename image_scanner_type,
        typename overlap_tester_type = test_box_overlap
        >
    class structural_object_detection_trainer : noncopyable
    {
        /*!
            REQUIREMENTS ON image_scanner_type
                image_scanner_type must be an implementation of 
                dlib/image_processing/scan_image_pyramid_abstract.h

            REQUIREMENTS ON overlap_tester_type
                overlap_tester_type must be an implementation of the test_box_overlap
                object defined in dlib/image_processing/box_overlap_testing_abstract.h.

            WHAT THIS OBJECT REPRESENTS
                This object is a tool for learning to detect objects in images based on a 
                set of labeled images. The training procedure produces an object_detector 
                which can be used to predict the locations of objects in new images.

                Note that this is just a convenience wrapper around the structural_svm_object_detection_problem 
                to make it look similar to all the other trainers in dlib.  
        !*/

    public:
        typedef double scalar_type;
        typedef default_memory_manager mem_manager_type;
        typedef object_detector<image_scanner_type,overlap_tester_type> trained_function_type;


        explicit structural_object_detection_trainer (
            const image_scanner_type& scanner
        );
        /*!
            requires
                - scanner.get_num_detection_templates() > 0
            ensures
                - #get_c() == 1
                - this object isn't verbose
                - #get_epsilon() == 0.3
                - #get_num_threads() == 2
                - #get_max_cache_size() == 40
                - #get_match_eps() == 0.5
                - #get_loss_per_missed_target() == 1
                - #get_loss_per_false_alarm() == 1
                - This object will attempt to learn a model for the given
                  scanner object when train() is called.
                - if (overlap_tester_type == test_box_overlap) then
                    - #auto_set_overlap_tester() == true
                - else
                    - #auto_set_overlap_tester() == false
        !*/

        bool auto_set_overlap_tester (
        ) const;
        /*!
            ensures
                - if (this object will automatically determine an appropriate 
                  state for the overlap tester used for non-max suppression.) then
                    - returns true
                    - In this case, it is determined using the find_tight_overlap_tester() 
                      routine based on the truth_rects given to the 
                      structural_object_detection_trainer::train() method.  
                - else
                    - returns false
        !*/

        void set_overlap_tester (
            const overlap_tester_type& tester
        );
        /*!
            ensures
                - #get_overlap_tester() == tester
                - #auto_set_overlap_tester() == false
        !*/

        overlap_tester_type get_overlap_tester (
        ) const;
        /*!
            requires
                - auto_set_overlap_tester() == false
            ensures
                - returns the overlap tester object which will be used to perform non-max suppression.
                  In particular, this function returns the overlap tester which will populate the
                  object_detector returned by train().
        !*/

        void set_num_threads (
            unsigned long num
        );
        /*!
            ensures
                - #get_num_threads() == num
        !*/

        unsigned long get_num_threads (
        ) const;
        /*!
            ensures
                - returns the number of threads used during training.  You should 
                  usually set this equal to the number of processing cores on your
                  machine.
        !*/

        void set_epsilon (
            scalar_type eps
        );
        /*!
            requires
                - eps > 0
            ensures
                - #get_epsilon() == eps
        !*/

        const scalar_type get_epsilon (
        ) const;
        /*!
            ensures
                - returns the error epsilon that determines when training should stop.
                  Smaller values may result in a more accurate solution but take longer 
                  to train.  You can think of this epsilon value as saying "solve the 
                  optimization problem until the average loss per sample is within epsilon 
                  of its optimal value".
        !*/

        void set_max_cache_size (
            unsigned long max_size
        );
        /*!
            ensures
                - #get_max_cache_size() == max_size
        !*/

        unsigned long get_max_cache_size (
        ) const;
        /*!
            ensures
                - During training, this object basically runs the object detector on 
                  each image, over and over.  To speed this up, it is possible to cache
                  the results of these detector invocations.  This function returns the 
                  number of cache elements per training sample kept in the cache.  Note 
                  that a value of 0 means caching is not used at all.  Note also that 
                  each cache element takes up about sizeof(double)*scanner.get_num_dimensions()
                  memory (where scanner is the scanner given to this object's constructor).
        !*/

        void be_verbose (
        );
        /*!
            ensures
                - This object will print status messages to standard out so that a 
                  user can observe the progress of the algorithm.
        !*/

        void be_quiet (
        );
        /*!
            ensures
                - this object will not print anything to standard out
        !*/

        void set_oca (
            const oca& item
        );
        /*!
            ensures
                - #get_oca() == item 
        !*/

        const oca get_oca (
        ) const;
        /*!
            ensures
                - returns a copy of the optimizer used to solve the structural SVM problem.  
        !*/

        void set_c (
            scalar_type C
        );
        /*!
            requires
                - C > 0
            ensures
                - #get_c() = C
        !*/

        const scalar_type get_c (
        ) const;
        /*!
            ensures
                - returns the SVM regularization parameter.  It is the parameter 
                  that determines the trade-off between trying to fit the training 
                  data (i.e. minimize the loss) or allowing more errors but hopefully 
                  improving the generalization of the resulting detector.  Larger 
                  values encourage exact fitting while smaller values of C may encourage 
                  better generalization. 
        !*/

        void set_match_eps (
            double eps
        );
        /*!
            requires
                - 0 < eps < 1
            ensures
                - #get_match_eps() == eps
        !*/

        double get_match_eps (
        ) const;
        /*!
            ensures
                - returns the amount of alignment necessary for a detection to be considered
                  as matching with a ground truth rectangle.  If it doesn't match then
                  it is considered to be a false alarm.  To define this precisely, let
                  A and B be two rectangles, then A and B match if and only if:
                    A.intersect(B).area()/(A+B).area() > get_match_eps()
        !*/

        double get_loss_per_missed_target (
        ) const;
        /*!
            ensures
                - returns the amount of loss experienced for failing to detect one of the
                  targets.  If you care more about finding targets than having a low false
                  alarm rate then you can increase this value.
        !*/

        void set_loss_per_missed_target (
            double loss
        );
        /*!
            requires
                - loss > 0
            ensures
                - #get_loss_per_missed_target() == loss
        !*/

        double get_loss_per_false_alarm (
        ) const;
        /*!
            ensures
                - returns the amount of loss experienced for emitting a false alarm detection.
                  Or in other words, the loss for generating a detection that doesn't correspond 
                  to one of the truth rectangles.  If you care more about having a low false
                  alarm rate than finding all the targets then you can increase this value.
        !*/

        void set_loss_per_false_alarm (
            double loss
        );
        /*!
            requires
                - loss > 0
            ensures
                - #get_loss_per_false_alarm() == loss
        !*/

        template <
            typename image_array_type
            >
        const trained_function_type train (
            const image_array_type& images,
            const std::vector<std::vector<rectangle> >& truth_rects
        ) const;
        /*!
            requires
                - is_learning_problem(images, truth_rects) == true
                - it must be valid to pass images[0] into the image_scanner_type::load() method.
                  (also, image_array_type must be an implementation of dlib/array/array_kernel_abstract.h)
            ensures
                - Uses the structural_svm_object_detection_problem to train an object_detector 
                  on the given images and truth_rects.  
                - returns a function F with the following properties:
                    - F(new_image) == A prediction of what objects are present in new_image.  This
                      is a set of rectangles indicating their positions.
        !*/
    }; 

// ----------------------------------------------------------------------------------------

}

#endif // DLIB_STRUCTURAL_OBJECT_DETECTION_TRAiNER_H_ABSTRACT__


