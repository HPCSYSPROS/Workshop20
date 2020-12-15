#include <cassert>
#include <cstdio>
#include <cstring>
#include <iostream>
#include <list>
#include <string>
#include <vector>

#define H5_USE_16_API 1
#include <cctk_Config.h>

#include <hdf5.h>

using namespace std;

// Global constants
int const dim = 3;

// Files and datasets
int const num_input_files = 0;
char const *const input_file_pattern = "interptoarray::parrays3d.file_%d.h5";
int const iteration_divisor = 1024;
char const *const output_file_name = "parrays3d-float.h5";
char const *const output_dataset_name = "phi";
hsize_t const output_dims[] = {0, 1024, 1024, 1024}; // [t,z,y,x]
hsize_t const output_maxdims[] = {H5S_UNLIMITED, 1024, 1024, 1024};
typedef float output_t;
hid_t const output_type = H5T_NATIVE_FLOAT;

// Technical details
hsize_t const chunk_size[] = {1, 128, 128, 128}; // 8 MB for float
int const compression_level = 1;
bool const write_checksum = true;

// Check a return value
#define check(_expr)                                                           \
  do {                                                                         \
    bool const _val = (_expr);                                                 \
    assert(_val);                                                              \
  } while (0)

static herr_t add_name(hid_t group, const char *name, H5L_info_t const *info,
                       void *op_data);

int main(int argc, char **argv) {
  cout << "hdf5_recombiner" << endl
       << "Copyright 2009 Erik Schnetter <schnetter@cct.lsu.edu>" << endl
       << endl;

  cout << "Opening output file \"" << output_file_name << "\"" << endl;

  htri_t is_hdf5;
  H5E_BEGIN_TRY { is_hdf5 = H5Fis_hdf5(output_file_name); }
  H5E_END_TRY;
  bool const file_exists = is_hdf5 > 0;
  hid_t const output_file =
      file_exists
          ? H5Fopen(output_file_name, H5F_ACC_RDWR, H5P_DEFAULT)
          : H5Fcreate(output_file_name, H5F_ACC_EXCL, H5P_DEFAULT, H5P_DEFAULT);

  cout << "Opening output dataset \"" << output_dataset_name << "\"" << endl;

  hid_t const output_datatype = output_type;

  hid_t const output_dataspace =
      H5Screate_simple(dim + 1, output_dims, output_maxdims);
  assert(output_dataspace >= 0);

  hid_t const output_properties = H5Pcreate(H5P_DATASET_CREATE);
  assert(output_properties >= 0);
  check(not H5Pset_chunk(output_properties, dim + 1, chunk_size));
  if (compression_level > 0) {
    check(not H5Pset_shuffle(output_properties));
    check(not H5Pset_deflate(output_properties, compression_level));
  }
  if (write_checksum) {
    check(not H5Pset_fletcher32(output_properties));
  }

  hid_t const output_dataset =
      file_exists ? H5Dopen(output_file, output_dataset_name)
                  : H5Dcreate(output_file, output_dataset_name, output_datatype,
                              output_dataspace, output_properties);
  assert(output_dataset >= 0);

  list<string> input_file_names;
  for (int input_file_num = 0; input_file_num < num_input_files;
       ++input_file_num) {
    char input_file_name[10000];
    snprintf(input_file_name, sizeof input_file_name, input_file_pattern,
             input_file_num);
    input_file_names.push_back(input_file_name);
  }
  for (int n = 1; n < argc; ++n) {
    input_file_names.push_back(argv[n]);
  }

  for (list<string>::const_iterator iinput_file_name = input_file_names.begin();
       iinput_file_name != input_file_names.end(); ++iinput_file_name) {
    string const &input_file_name = *iinput_file_name;
    cout << "Opening input file \"" << input_file_name << "\"" << endl;
    hid_t const input_file =
        H5Fopen(input_file_name.c_str(), H5F_ACC_RDWR, H5P_DEFAULT);
    assert(input_file >= 0);

    typedef list<string> names_t;
    names_t names;
    hsize_t idx = 0;
    check(not H5Literate(input_file, H5_INDEX_NAME, H5_ITER_NATIVE, &idx,
                         add_name, &names));
    for (names_t::const_iterator iname = names.begin(); iname != names.end();
         ++iname) {
      char const *const input_dataset_name = iname->c_str();
      cout << "Reading input dataset \"" << input_dataset_name << "\"" << endl;

      char varname[10000];
      int it, tl, c;
      int const icnt = sscanf(input_dataset_name, "%s it=%d tl=%d c=%d",
                              varname, &it, &tl, &c);
      assert(icnt == 4);
      assert(it % iteration_divisor == 0);
      int const iteration = it / iteration_divisor;

      hid_t const input_dataset = H5Dopen(input_file, input_dataset_name);
      assert(input_dataset >= 0);

      int iorigin[dim];
      hid_t const iorigin_attr = H5Aopen(input_dataset, "iorigin", H5P_DEFAULT);
      assert(iorigin_attr >= 0);
      check(not H5Aread(iorigin_attr, H5T_NATIVE_INT, iorigin));
      check(not H5Aclose(iorigin_attr));

      hid_t const input_dataspace = H5Dget_space(input_dataset);
      assert(input_dataspace >= 0);
      assert(H5Sis_simple(input_dataspace) > 0);
      assert(H5Sget_simple_extent_ndims(input_dataspace) == dim);
      hsize_t input_dims[dim];
      check(H5Sget_simple_extent_dims(input_dataspace, input_dims, NULL) ==
            dim);

      hid_t const input_memory_dataspace =
          H5Screate_simple(dim, input_dims, NULL);
      assert(input_memory_dataspace);
      hssize_t const input_memory_npoints =
          H5Sget_simple_extent_npoints(input_memory_dataspace);
      vector<output_t> data(input_memory_npoints);

      check(not H5Dread(input_dataset, output_type, input_memory_dataspace,
                        input_dataspace, H5P_DEFAULT, &data.front()));

      check(not H5Sclose(input_memory_dataspace));
      check(not H5Sclose(input_dataspace));
      check(not H5Dclose(input_dataset));

      cout << "Writing output dataset" << endl;

      hsize_t output_extent[dim + 1];
      output_extent[0] = iteration + 1;
      for (int d = 0; d < dim; ++d) {
        output_extent[d + 1] = output_dims[d + 1];
      }
      check(not H5Dextend(output_dataset, output_extent));

      hid_t const output_dataspace =
          H5Screate_simple(dim + 1, output_extent, output_maxdims);
      assert(output_dataspace >= 0);

      hsize_t output_offset[dim + 1];
      hsize_t output_dims[dim + 1];
      output_offset[0] = iteration;
      output_dims[0] = 1;
      for (int d = 0; d < dim; ++d) {
        output_offset[d + 1] = iorigin[dim - 1 - d];
        output_dims[d + 1] = input_dims[d];
      }
      check(not H5Sselect_hyperslab(output_dataspace, H5S_SELECT_SET,
                                    output_offset, NULL, output_dims, NULL));
      cout << "   extent [" << output_offset[3] << "," << output_offset[2]
           << "," << output_offset[1] << "," << output_offset[0] << "] - ["
           << output_offset[3] + output_dims[3] << ","
           << output_offset[2] + output_dims[2] << ","
           << output_offset[1] + output_dims[1] << ","
           << output_offset[0] + output_dims[0] << "] "
           << "(" << data.size() << " points)" << endl;

      hid_t const output_memory_dataspace =
          H5Screate_simple(dim + 1, output_dims, NULL);
      assert(output_memory_dataspace);
      hssize_t const output_memory_npoints =
          H5Sget_simple_extent_npoints(output_memory_dataspace);
      assert(output_memory_npoints == input_memory_npoints);

      check(not H5Dwrite(output_dataset, output_type, output_memory_dataspace,
                         output_dataspace, H5P_DEFAULT, &data.front()));

      H5Sclose(output_memory_dataspace);
      H5Sclose(output_dataspace);

    } // for names

    check(not H5Fclose(input_file));

  } // for input_file_num

  check(not H5Dclose(output_dataset));
  check(not H5Pclose(output_properties));
  check(not H5Sclose(output_dataspace));
  check(not H5Fclose(output_file));

  cout << "Done." << endl;
  return 0;
}

static herr_t add_name(hid_t const group, const char *const name,
                       H5L_info_t const *const info, void *const op_data) {
  typedef list<string> names_t;
  names_t &names = *static_cast<names_t *>(op_data);
  if (strcmp(name, "Parameters and Global Attributes") != 0) {
    names.push_back(name);
  }
  return 0;
}
