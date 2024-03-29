/*
 Copyright (C) 2019 X. Andrade, 2021 M. Lueders

 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation; either version 2, or (at your option)
 any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program; if not, write to the Free Software
 Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA
 02110-1301, USA.

*/

#include <config.h>

#include <cassert>
#include <string>
#include <unordered_map>

typedef std::unordered_map<std::string, void *> map_type;

extern "C" void sphash_map_init(map_type **map) { *map = new map_type; }

extern "C" void sphash_map_end(map_type **map) { delete *map; }

extern "C" void sphash_map_insert(map_type *map, const char *key, void *ptr) {
  std::string my_key(key);
  assert(map->count(my_key) == 0);
  (*map)[my_key] = ptr;
}

extern "C" void sphash_map_lookup(const map_type *map, const char *key,
                                  int *found, void **ptr) {
  std::string my_key(key);
  auto it = (map)->find(my_key);

  if (it == (map)->end()) {
    *found = 0;
  } else {
    *found = 1;
    *ptr = it->second;
  }
}

extern "C" void sphash_iterator_low_start(map_type::const_iterator *iterator,
                                          map_type::const_iterator *end,
                                          const map_type *map) {
  *iterator = map->begin();
  *end = map->end();
}

extern "C" void sphash_iterator_low_has_next(map_type::const_iterator iterator,
                                             map_type::const_iterator end,
                                             int *result) {
  *result = (iterator != end);
}

extern "C" void sphash_iterator_low_get(map_type::const_iterator *iterator,
                                        void **ptr) {
  *ptr = (*iterator)->second;
  (*iterator)++;
}
