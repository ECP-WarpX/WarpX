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

#include <AMReX_Vector.H>

#include <map>
#include <string>
#include <type_traits>


namespace ablastr::utils
{
    /** This is a std::map where operator[] throws if the key is missing.
     *
     * This also supports T* to T const* conversion for T.
     */
    template< typename Key, typename T >
    struct ConstMap : public std::map< Key, T >
    {
        // default constructors
        ConstMap () = default;
        ConstMap (ConstMap const &) = default;
        ConstMap (ConstMap &&) = default;
        ConstMap& operator=(ConstMap const & other) = default;
        ConstMap& operator=(ConstMap && other) noexcept = default;

        // conversion to const value type
        template<
            typename mutable_value = std::remove_const_t<std::remove_pointer_t<T>>,
            std::enable_if_t<
                !std::is_same<mutable_value, T>{},
                int
            > = 0
        >
        ConstMap (ConstMap<Key, mutable_value*> const & other) {
            for (auto const & kv : other) {
                this->insert({kv.first, kv.second});
            }
        }

        // redefine operators
        T operator[]( const Key& key )
        {
            return this->at(key);
        }

        T operator[]( Key&& key )
        {
            return this->at(key);
        }

        T operator[]( const Key& key ) const
        {
            return this->at(key);
        }

        T operator[]( Key&& key ) const
        {
            return this->at(key);
        }
    };

    /** Advanced vector that supports T to T const conversion
     */
    template< typename T >
    struct MultiLevelVect : public amrex::Vector< T >
    {
        // default constructors
        MultiLevelVect () = default;
        MultiLevelVect (MultiLevelVect const &) = default;
        MultiLevelVect (MultiLevelVect &&) = default;
        MultiLevelVect& operator=(MultiLevelVect const & other) = default;
        MultiLevelVect& operator=(MultiLevelVect && other) noexcept = default;

        // conversion to const value type
        template<
            typename mutable_value = std::remove_const_t<T>,
            std::enable_if_t<
                !std::is_same<mutable_value, T>{},
                int
            > = 0
        >
        MultiLevelVect (MultiLevelVect<mutable_value> const & other) {
            for (auto const & v : other) {
                this->push_back(v);
            }
        }
    };

}  // ablastr::utils

#endif // ABLASTR_CONST_MAP_H_
