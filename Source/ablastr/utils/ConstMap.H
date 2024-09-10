/* Copyright 2024 The ABLASTR Community
 *
 * This file is part of ABLASTR.
 *
 * License: BSD-3-Clause-LBNL
 * Authors: Axel Huebl
 */

#ifndef ABLASTR_CONST_MAP_H_
#define ABLASTR_CONST_MAP_H_

#include "TextMsg.H"

#include <map>
#include <string>


namespace ablastr::utils
{
    /** This is a std::map where operator[] throws if the key is missing
     */
    template< typename Key, typename T >
    struct ConstMap : public std::map< Key, T >
    {
        T& operator[]( const Key& key )
        {
            return this->at(key);
        }

        T& operator[]( Key&& key )
        {
            return this->at(key);
        }

        T& operator[]( const Key& key ) const
        {
            return this->at(key);
        }

        T& operator[]( Key&& key ) const
        {
            return this->at(key);
        }
    };

}  // ablastr::utils

#endif // ABLASTR_CONST_MAP_H_