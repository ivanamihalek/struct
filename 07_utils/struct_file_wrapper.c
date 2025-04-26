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

/* ======================== WRAPPER IMPLEMENTATION ======================== */
#include <errno.h>
#include <stddef.h>  //otherwise clion complains it cannot finf NULL
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <zlib.h>

#include "struct.h"
#include "struct_utils.h"


// Error reporting
const char* file_error(const FileWrapper* fw) { return fw->error_msg; }
int file_errno(const FileWrapper* fw) { return fw->last_error; }
void file_clearerr(FileWrapper* fw) { fw->last_error = 0; fw->error_msg = NULL; }

FileWrapper* file_open(const char* path, const char* mode) {
    FileWrapper* fw = calloc(1, sizeof(FileWrapper));
    if (!fw) return NULL;

    const int gzipped = is_gzipped(path);
    if (gzipped == ERR_NO_FILE_OR_CHAIN) {
        fprintf(stderr, "File \"%s\" does not exist\n", path);
        exit(ERR_NO_FILE_OR_CHAIN);
    }
    fw->comp_type = gzipped ? COMP_GZIP : COMP_NONE;
    file_clearerr(fw);

    if (fw->comp_type == COMP_NONE) {
        fw->handle.file = fopen(path, mode);
        if (!fw->handle.file) {
            fprintf(stderr, "Errno %d opening plain file \"%s\" \n", errno, path);
            exit(errno);  // a macro defined in errno.h
        }
    } else {
        fw->handle.gz = gzopen(path, mode);
        if (!fw->handle.gz) {
            return NULL;
            fprintf(stderr, "Errno %d opening gzipped file \"%s\" \n", Z_ERRNO, path);
            exit(Z_ERRNO);  // a macro defined in errno.h
        }
    }
    return fw;
}

ssize_t file_read(void* buf, size_t size, FileWrapper* fw) {
    file_clearerr(fw);

    if (fw->comp_type == COMP_NONE) {
        size_t read = fread(buf, 1, size, fw->handle.file);
        if (ferror(fw->handle.file)) {
            fw->last_error = errno;
            fw->error_msg = strerror(errno);
            return -1;
        }
        return (ssize_t)read;
    }

    int zread = gzread(fw->handle.gz, buf, size);
    if (zread < 0) {
        int zerr;
        fw->error_msg = gzerror(fw->handle.gz, &zerr);
        fw->last_error = zerr;
        return -1;
    }
    return zread;
}

int file_seek(FileWrapper* fw, long offset, int whence) {
    file_clearerr(fw);

    if (fw->comp_type == COMP_NONE) {
        int ret = fseek(fw->handle.file, offset, whence);
        if (ret != 0) {
            fw->last_error = errno;
            fw->error_msg = strerror(errno);
        }
        return ret;
    }

    z_off_t ret = gzseek(fw->handle.gz, offset, whence);
    if (ret == -1) {
        int zerr;
        fw->error_msg = gzerror(fw->handle.gz, &zerr);
        fw->last_error = zerr;
        return -1;
    }
    return 0;
}

int file_rewind(FileWrapper* fw) {
    file_clearerr(fw);

    if (fw->comp_type == COMP_NONE) {
        rewind(fw->handle.file);
        if (ferror(fw->handle.file)) {
            fw->last_error = errno;
            fw->error_msg = strerror(errno);
            return -1;
        }
        return 0;
    }

    int ret = gzrewind(fw->handle.gz);
    if (ret != 0) {
        int zerr;
        fw->error_msg = gzerror(fw->handle.gz, &zerr);
        fw->last_error = zerr;
        return -1;
    }
    return 0;
}

int file_feof(FileWrapper* fw) {
    if (fw->comp_type == COMP_NONE) {
        return feof(fw->handle.file);
    }
    return gzeof(fw->handle.gz);
}

char* file_gets(char* buf, int size, FileWrapper* fw) {
    file_clearerr(fw);

    if (fw->comp_type == COMP_NONE) {
        char* ret = fgets(buf, size, fw->handle.file);
        if (!ret) {
            if (ferror(fw->handle.file)) {
                fw->last_error = errno;
                fw->error_msg = strerror(errno);
            }
            return NULL;
        }
        return buf;
    }

    char* ret = gzgets(fw->handle.gz, buf, size);
    if (!ret) {
        if (!gzeof(fw->handle.gz)) {  // Only set error if not EOF
            int zerr;
            fw->error_msg = gzerror(fw->handle.gz, &zerr);
            fw->last_error = zerr;
        }
        return NULL;
    }

    // Strip trailing newline if present (matches fgets behavior)
    size_t len = strlen(buf);
    if (len > 0 && buf[len-1] == '\n') {
        buf[len-1] = '\0';
    }

    return buf;
}

void file_close(FileWrapper* fw) {
    if (fw) {
        if (fw->comp_type == COMP_NONE) {
            if (fw->handle.file) fclose(fw->handle.file);
        } else {
            if (fw->handle.gz) gzclose(fw->handle.gz);
        }
        free(fw);
    }
}
