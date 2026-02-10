#include "output_file_manager.h"
#include "utilities.h"
#include "definitions.h"

OutputFileManager::OutputFileManager()
  : vm_file (nullptr)
  , str_str_file (nullptr)
  , err_file (nullptr)
  , conv_file (nullptr)
{
  // open_files if my_rank == 0 form EVPFFT constructor
  //open_files();
}

void OutputFileManager::open_files()
{
#ifndef ABSOLUTE_NO_OUTPUT
    vm_file = std::shared_ptr <FILE> (fopen("vm.out", "w"), &fclose);
    check_that_file_is_open(vm_file.get(), "vm.out");

    str_str_file = std::shared_ptr <FILE> (fopen("str_str.out", "w"), &fclose);
    check_that_file_is_open(str_str_file.get(), "str_str.out");

    err_file = std::shared_ptr <FILE> (fopen("err.out", "w"), &fclose);
    check_that_file_is_open(err_file.get(), "err.out");

    conv_file = std::shared_ptr <FILE> (fopen("conv.out", "w"), &fclose);
    check_that_file_is_open(conv_file.get(), "conv.out");


    // write heading text
    fprintf(vm_file.get(), "EVM,EVMP,DVM,DVMP,SVM,SVM1\n");

    fprintf(str_str_file.get(),
            "E11,E22,E33,E23,E13,E12,"
            "S11,S22,S33,S23,S13,S12,"
            "EP11,EP22,EP33,EP23,EP13,EP12,"
            "E_dot11,E_dot22,E_dot33,E_dot23,E_dot13,E_dot12,"
            "E_dotP11,E_dotP22,E_dotP33,E_dotP23,E_dotP13,E_dotP12,"
            "EVM,EPVM,DVM,DPVM,SVM\n");

    fprintf(err_file.get(), "STEP,IT,ERRE,ERRS,SVM,AVG_NR_IT\n");
#endif
}
