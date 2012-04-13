// Copyright (C) 2011  Davis E. King (davis@dlib.net)
// License: Boost Software License   See LICENSE.txt for the full license.

#include <sstream>
#include <string>
#include <cstdlib>
#include <ctime>
#include <dlib/hash.h>
#include <dlib/matrix.h>
#include <dlib/byte_orderer.h>

#include "tester.h"

namespace  
{
    using namespace test;
    using namespace dlib;
    using namespace std;

    logger dlog("test.hash");


    template <typename T>
    void to_little (
        std::vector<T>& item
    )
    {
        byte_orderer bo;
        for (unsigned long i = 0; i < item.size(); ++i)
            bo.host_to_little(item[i]);
    }


    template <typename T>
    void to_little (
        matrix<T>& item
    )
    {
        byte_orderer bo;
        for (long r = 0; r < item.nr(); ++r)
        {
            for (long c = 0; c < item.nc(); ++c)
            {
                bo.host_to_little(item(r,c));
            }
        }
    }

    // Run the official test for MurmurHash3
    void murmur_hash_test()
    {
        uint8 key[256];
        uint32 hashes[256];
        uint32 final = 0;

        memset(key,0,sizeof(key));
        memset(hashes,0,sizeof(hashes));

        // Hash keys of the form {0}, {0,1}, {0,1,2}... up to N=255,using 256-N as
        // the seed.
        for(int i = 0; i < 256; i++)
        {
            key[i] = (uint8)i;

            hashes[i] = murmur_hash3(key,i,256-i);
        }

        byte_orderer bo;
        bo.host_to_little(hashes);
        final = murmur_hash3(hashes,sizeof(hashes),0);

        // using ostringstream to avoid compiler error in visual studio 2005
        ostringstream sout;
        sout << hex << final;
        dlog << LINFO << "final: "<< sout.str();
        DLIB_TEST(final == 0xB0F57EE3);
    }

    void murmur_hash_128_test()
    {
        uint8 key[256];
        uint64 hashes[256*2];
        uint32 final = 0;

        memset(key,0,sizeof(key));
        memset(hashes,0,sizeof(hashes));

        // Hash keys of the form {0}, {0,1}, {0,1,2}... up to N=255,using 256-N as
        // the seed.
        for(int i = 0; i < 256; i++)
        {
            key[i] = (uint8)i;

            const std::pair<uint64,uint64> temp = murmur_hash3_128bit(key,i,256-i);
            hashes[2*i]   = temp.first;
            hashes[2*i+1] = temp.second;
        }

        byte_orderer bo;
        bo.host_to_little(hashes);
        final = static_cast<uint32>(murmur_hash3_128bit(hashes,sizeof(hashes),0).first);

        // using ostringstream to avoid compiler error in visual studio 2005
        ostringstream sout;
        sout << hex << final;
        dlog << LINFO << "final 64: "<< sout.str();
        DLIB_TEST(final == 0x6384BA69);
    }

    class test_hash : public tester
    {
    public:
        test_hash (
        ) :
            tester ("test_hash",
                    "Runs tests on the hash routines.")
        {}

        void perform_test (
        )
        {
            print_spinner();

            murmur_hash_test();
            murmur_hash_128_test();

            std::string str1 = "some random string";
            matrix<unsigned char> mat(2,2);

            mat = 1,2,3,4;

            matrix<uint64> mat2(2,3);

            mat2 = 1,2,3,4,5,6;

            to_little(mat2);

            std::vector<unsigned char> v(4);
            v[0] = 'c';
            v[1] = 'a';
            v[2] = 't';
            v[3] = '!';

            std::vector<uint16> v2(4);
            v[0] = 'c';
            v[1] = 'a';
            v[2] = 't';
            v[3] = '!';
            to_little(v2);

            std::map<unsigned char, unsigned char> m;
            m['c'] = 'C';
            m['a'] = 'A';
            m['t'] = 'T';

            dlog << LINFO << "hash(str1): "<< dlib::hash(str1);
            dlog << LINFO << "hash(v):    "<< dlib::hash(v);
            dlog << LINFO << "hash(v2):   "<< dlib::hash(v2);
            dlog << LINFO << "hash(m):    "<< dlib::hash(m);
            dlog << LINFO << "hash(mat):  "<< dlib::hash(mat);
            dlog << LINFO << "hash(mat2): "<< dlib::hash(mat2);

            DLIB_TEST(dlib::hash(str1) == 0x3ffe6bf6);
            DLIB_TEST(dlib::hash(v)    == 0xf1af2ca6);
            DLIB_TEST(dlib::hash(v2)   == 0x63852afc);
            DLIB_TEST(dlib::hash(m)    == 0xaacc3f6f);
            DLIB_TEST(dlib::hash(mat)  == 0x3e349da5);
            DLIB_TEST(dlib::hash(mat2) == 0x3a95dc52);
            DLIB_TEST(murmur_hash3(&str1[0], str1.size(), 0) == 0x3ffe6bf6);

            dlog << LINFO << "hash(str1,1): "<< dlib::hash(str1,1);
            dlog << LINFO << "hash(v,3):    "<< dlib::hash(v,3);
            dlog << LINFO << "hash(v2,3):   "<< dlib::hash(v2,3);
            dlog << LINFO << "hash(m,4):    "<< dlib::hash(m,4);
            dlog << LINFO << "hash(mat,5):  "<< dlib::hash(mat,5);
            dlog << LINFO << "hash(mat2,6): "<< dlib::hash(mat2,6);

            DLIB_TEST(dlib::hash(str1,1) == 0xb17cea93);
            DLIB_TEST(dlib::hash(v,3)    == 0x7ec9284c);
            DLIB_TEST(dlib::hash(v2,3)   == 0xb2ce147f);
            DLIB_TEST(dlib::hash(m,4)    == 0xfa5e7ac2);
            DLIB_TEST(dlib::hash(mat,5)  == 0x8de27259);
            DLIB_TEST(dlib::hash(mat2,6) == 0xb8aa7714);
            DLIB_TEST(murmur_hash3(&str1[0], str1.size(), 1) == 0xb17cea93);

        }
    } a;



}



