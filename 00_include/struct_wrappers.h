/*
This source code is part of deconSTRUCT,
protein structure database search and backbone alignment application.
Copyright (C) 2008-2025 Ivana Mihalek.

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or,
at your option, any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program. If not, see <http://www.gnu.org/licenses/>.

Contact: ivana.mihalek@gmail.com.
*/

# ifndef STRUCT_WRAPPERS_H
# define STRUCT_WRAPPERS_H

# include <stdio.h>
# include <string.h>
# include <zlib.h>

typedef enum {
    NONE,
    GZIP
} COMPRESSION_TYPE;

// Define a function pointer type for reading
typedef int (*read_func_t)(char *buffer, int size, void* file);
// Define function pointer types for seeking and rewinding
typedef int (*seek_func_t)(void *file, long offset, int whence);
typedef void (*rewind_func_t)(void *file);

int read_wrapper(char* buffer, int size, void* file, COMPRESSION_TYPE compression_type);
int seek_wrapper(void *file, long offset, int whence, COMPRESSION_TYPE compression_type);
void rewind_wrapper(void *file, COMPRESSION_TYPE compression_type);

#endif //STRUCT_WRAPPERS_H
