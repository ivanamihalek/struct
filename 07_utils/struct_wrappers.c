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

# include <stddef.h>
# include <stdlib.h>
# include "struct_wrappers.h"

// Wrapper for fgets (regular files)
int read_wrapper(char* buffer, int size, void* file, COMPRESSION_TYPE compression_type) {
     switch (compression_type){
        case NONE:
            return fgets(buffer, size, (FILE *)file) != NULL ? strlen(buffer) : 0;
        case GZIP:
             return gzread((gzFile)file, buffer, size);
        default:
             fprintf(stderr, "Unrecognized compression type.")
             exit(1)
    }
}

// Wrapper for fgets (regular files)
int seek_wrapper(void *file, long offset, int whence, COMPRESSION_TYPE compression_type) {
     switch (compression_type){
        case NONE:
            return fseek((FILE *)file, offset, whence);
        case GZIP:
            return gzseek((gzFile)file, offset, whence);
        default:
             fprintf(stderr, "Unrecognized compression type.")
             exit(1)
    }
}

void rewind_wrapper(void *file, COMPRESSION_TYPE compression_type) {
     switch (compression_type){
        case NONE:
            return rewind((FILE *)file);
        case GZIP:
            return gzrewind((gzFile)file);
        default:
             fprintf(stderr, "Unrecognized compression type.")
             exit(1)
    }
}

