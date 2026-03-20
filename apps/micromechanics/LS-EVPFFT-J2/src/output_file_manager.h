#pragma once

#include <stdio.h>
#include <string>
#include <memory>



struct OutputFileManager
{

  std::shared_ptr <FILE> vm_file;
  std::shared_ptr <FILE> str_str_file;
  std::shared_ptr <FILE> err_file;
  std::shared_ptr <FILE> conv_file;

  OutputFileManager();
  void open_files();
};
