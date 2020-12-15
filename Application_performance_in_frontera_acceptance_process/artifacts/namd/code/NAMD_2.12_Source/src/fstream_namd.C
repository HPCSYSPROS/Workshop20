
#include "fstream_namd.h"
#include "common.h"

void ofstream_namd::open(const char *_fname, const char *_ext) {
  if ( fd ) NAMD_bug("ofstream_namd::open() called when file is already open");
  NAMD_backup_file(_fname, _ext);
  open(_fname);
}

void ofstream_namd::open(const char *_fname, std::ios_base::openmode _mode) {
  if ( fd ) NAMD_bug("ofstream_namd::open() called when file is already open");
  fname = _fname;
  fd = NAMD_open_text(_fname, _mode & std::ios_base::app);
}

ofstream_namd& ofstream_namd::flush() {
  if ( ! fd ) NAMD_bug("ofstream_namd::flush() called when file is not open");
  const std::string text = str();
  NAMD_write(fd, text.c_str(), text.size(), fname.c_str());
  str("");
  return *this;
}

void ofstream_namd::close() {
  if ( ! fd ) NAMD_bug("ofstream_namd::close() called when file is not open");
  flush();
  NAMD_close(fd, fname.c_str());
  fd = 0;
}

