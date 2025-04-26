//
// Created by ivana on 4/26/25.
//

#ifndef STRUCT_FILE_WRAPPER_H
#define STRUCT_FILE_WRAPPER_H

#include <zlib.h>

typedef enum { COMP_NONE, COMP_GZIP } CompressionType;

typedef struct {
    CompressionType comp_type;
    union {
        FILE* file;
        gzFile gz;
    } handle;
    int last_error;
    const char* error_msg;
} FileWrapper;

FileWrapper* file_open(const char* path, const char* mode);
ssize_t file_read(void* buf, size_t size, FileWrapper* fw);
int file_feof(FileWrapper* fw);
int file_seek(FileWrapper* fw, long offset, int whence);
int file_rewind(FileWrapper* fw);
void file_close(FileWrapper* fw);

char* file_gets(char* buf, int size, FileWrapper* fw);

#endif //STRUCT_FILE_WRAPPER_H
