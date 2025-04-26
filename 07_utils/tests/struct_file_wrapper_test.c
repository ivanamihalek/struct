//
// Created by ivana on 3/8/25.
//

/*
# Compile with:
gcc -Wall -o file_wrapper_test file_wrapper.c -lz -lcunit -DUNIT_TEST

# Run tests:
./file_wrapper_test

CUnit - A unit testing framework for C - Version 2.1-3
     ......
     Run Summary:    Type  Total    Ran Passed Failed Inactive
              suites      1      1    n/a      0        0
              tests       4      4      4      0        0
             asserts     19     19     19      0      n/a

*/

#ifdef UNIT_TEST
#include <CUnit/CUnit.h>
#include <CUnit/Basic.h>

#define TEST_TEXT "The quick brown fox jumps over the lazy dog\n"

static void create_test_files() {
    // Create uncompressed test file
    FILE* f = fopen("test.txt", "w");
    fputs(TEST_TEXT, f);
    fclose(f);

    // Create compressed test file
    gzFile gz = gzopen("test.gz", "wb");
    gzwrite(gz, TEST_TEXT, sizeof(TEST_TEXT)-1); // Exclude null terminator
    gzclose(gz);
}

static void remove_test_files() {
    remove("test.txt");
    remove("test.gz");
}

void test_file_open_close(void) {
    // Test valid files
    FileWrapper* f1 = file_open("test.txt", "rb", COMP_NONE);
    CU_ASSERT_PTR_NOT_NULL(f1);
    CU_ASSERT_EQUAL(f1->type, COMP_NONE);
    file_close(f1);

    FileWrapper* f2 = file_open("test.gz", "rb", COMP_GZIP);
    CU_ASSERT_PTR_NOT_NULL(f2);
    CU_ASSERT_EQUAL(f2->type, COMP_GZIP);
    file_close(f2);

    // Test invalid files
    FileWrapper* f3 = file_open("nonexistent.txt", "rb", COMP_NONE);
    CU_ASSERT_PTR_NULL(f3);

    FileWrapper* f4 = file_open("nonexistent.gz", "rb", COMP_GZIP);
    CU_ASSERT_PTR_NULL(f4);
}

void test_file_read(void) {
    char buffer[256];

    FileWrapper* f = file_open("test.txt", "rb", COMP_NONE);
    CU_ASSERT_PTR_NOT_NULL(f);

    ssize_t read = file_read(buffer, sizeof(buffer), f);
    CU_ASSERT_EQUAL(read, strlen(TEST_TEXT));
    CU_ASSERT_STRING_EQUAL(buffer, TEST_TEXT);
    file_close(f);

    FileWrapper* gzf = file_open("test.gz", "rb", COMP_GZIP);
    CU_ASSERT_PTR_NOT_NULL(gzf);

    read = file_read(buffer, sizeof(buffer), gzf);
    CU_ASSERT_EQUAL(read, strlen(TEST_TEXT));
    CU_ASSERT_STRING_EQUAL(buffer, TEST_TEXT);
    file_close(gzf);
}

void test_file_seek_rewind(void) {
    char buffer[256];

    // Test uncompressed file
    FileWrapper* f = file_open("test.txt", "rb", COMP_NONE);
    CU_ASSERT_PTR_NOT_NULL(f);

    CU_ASSERT_EQUAL(file_seek(f, 4, SEEK_SET), 0);
    CU_ASSERT_EQUAL(file_read(buffer, 5, f), 5);
    CU_ASSERT_STRING_EQUAL(buffer, "quick");

    CU_ASSERT_EQUAL(file_rewind(f), 0);
    CU_ASSERT_EQUAL(file_read(buffer, 3, f), 3);
    CU_ASSERT_STRING_EQUAL(buffer, "The");
    file_close(f);

    // Test compressed file
    FileWrapper* gzf = file_open("test.gz", "rb", COMP_GZIP);
    CU_ASSERT_PTR_NOT_NULL(gzf);

    CU_ASSERT_EQUAL(file_seek(gzf, 4, SEEK_SET), 0);
    CU_ASSERT_EQUAL(file_read(buffer, 5, gzf), 5);
    CU_ASSERT_STRING_EQUAL(buffer, "quick");

    CU_ASSERT_EQUAL(file_rewind(gzf), 0);
    CU_ASSERT_EQUAL(file_read(buffer, 3, gzf), 3);
    CU_ASSERT_STRING_EQUAL(buffer, "The");
    file_close(gzf);
}

void test_error_handling(void) {
    FileWrapper* f = file_open("test.txt", "rb", COMP_NONE);
    file_close(f);  // Close immediately

    // Test read on closed file
    char buffer[10];
    CU_ASSERT_EQUAL(file_read(buffer, 10, f), -1);
}

void test_file_gets(void) {
    char buffer[100];

    // Test uncompressed file
    FileWrapper* f = file_open("test.txt", "rb", COMP_NONE);
    CU_ASSERT_PTR_NOT_NULL(f);

    CU_ASSERT_PTR_NOT_NULL(file_gets(buffer, sizeof(buffer), f));
    CU_ASSERT_STRING_EQUAL(buffer, "The quick brown fox jumps over the lazy dog");
    file_close(f);

    // Test compressed file
    FileWrapper* gzf = file_open("test.gz", "rb", COMP_GZIP);
    CU_ASSERT_PTR_NOT_NULL(gzf);

    CU_ASSERT_PTR_NOT_NULL(file_gets(buffer, sizeof(buffer), gzf));
    CU_ASSERT_STRING_EQUAL(buffer, "The quick brown fox jumps over the lazy dog");
    file_close(gzf);

    // Test buffer limits
    FileWrapper* f2 = file_open("test.txt", "rb", COMP_NONE);
    CU_ASSERT_PTR_NOT_NULL(file_gets(buffer, 10, f2));
    CU_ASSERT_STRING_EQUAL(buffer, "The quic");
    file_close(f2);
}


void test_file_feof(void) {
    char buffer[256];

    // Test uncompressed file
    FileWrapper* f = file_open("test.txt", "rb", COMP_NONE);
    CU_ASSERT_PTR_NOT_NULL(f);

    // Read entire file
    while (file_read(buffer, sizeof(buffer), f) > 0);
    CU_ASSERT(file_feof(f));
    file_close(f);

    // Test compressed file
    FileWrapper* gzf = file_open("test.gz", "rb", COMP_GZIP);
    CU_ASSERT_PTR_NOT_NULL(gzf);

    // Read entire file
    while (file_read(buffer, sizeof(buffer), gzf) > 0);
    CU_ASSERT(file_feof(gzf));
    file_close(gzf);

    // Test mid-file reading
    FileWrapper* f2 = file_open("test.txt", "rb", COMP_NONE);
    CU_ASSERT_PTR_NOT_NULL(f2);
    CU_ASSERT_EQUAL(file_read(buffer, 10, f2), 10);
    CU_ASSERT(!file_feof(f2));
    file_close(f2);
}

int main() {
    CU_pSuite suite;

    CU_initialize_registry();
    suite = CU_add_suite("File Wrapper Tests", NULL, NULL);
    CU_add_test(suite, "Open/Close", test_file_open_close);
    CU_add_test(suite, "Read Operations", test_file_read);
    CU_add_test(suite, "Seek/Rewind", test_file_seek_rewind);
    CU_add_test(suite, "Error Handling", test_error_handling);
    CU_add_test(suite, "Fgets test, test_file_gets);
    CU_add_test(suite, "EOF Detection", test_file_feof);

    create_test_files();
    CU_basic_set_mode(CU_BRM_VERBOSE);
    CU_basic_run_tests();
    remove_test_files();

    CU_cleanup_registry();
    return CU_get_error();
}
#endif
