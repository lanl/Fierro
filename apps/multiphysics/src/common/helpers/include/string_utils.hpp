/**********************************************************************************************
� 2020. Triad National Security, LLC. All rights reserved.
This program was produced under U.S. Government contract 89233218CNA000001 for Los Alamos
National Laboratory (LANL), which is operated by Triad National Security, LLC for the U.S.
Department of Energy/National Nuclear Security Administration. All rights in the program are
reserved by Triad National Security, LLC, and the U.S. Department of Energy/National Nuclear
Security Administration. The Government is granted for itself and others acting on its behalf a
nonexclusive, paid-up, irrevocable worldwide license in this material to reproduce, prepare
derivative works, distribute copies to the public, perform publicly and display publicly, and
to permit others to do so.
This program is open source under the BSD-3 License.
Redistribution and use in source and binary forms, with or without modification, are permitted
provided that the following conditions are met:
1.  Redistributions of source code must retain the above copyright notice, this list of
conditions and the following disclaimer.
2.  Redistributions in binary form must reproduce the above copyright notice, this list of
conditions and the following disclaimer in the documentation and/or other materials
provided with the distribution.
3.  Neither the name of the copyright holder nor the names of its contributors may be used
to endorse or promote products derived from this software without specific prior
written permission.
THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR
CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
**********************************************************************************************/
#ifndef STRING_UTILITIES_H
#define STRING_UTILITIES_H

#include <string>
#include <vector>

extern const std::string WHITESPACE;

/////////////////////////////////////////////////////////////////////////////
///
/// \fn split
///
/// \brief Splits a string by a given delimiter
///
/// \param Input string
/// \param delimiter
///
/// \return Vector of split string values
///
/////////////////////////////////////////////////////////////////////////////
std::vector<std::string> split(std::string s, std::string delimiter);

/////////////////////////////////////////////////////////////////////////////
///
/// \fn exact_array_values
///
/// \brief Extracts array values between [ and ] (generally for yaml parsing)
///
/// \param Input string
/// \param delimiter
///
/// \return Vector of split string values
///
/////////////////////////////////////////////////////////////////////////////
std::vector<std::string> exact_array_values(std::string s, std::string delimiter);

/////////////////////////////////////////////////////////////////////////////
///
/// \fn exact_array_values
///
/// \brief Extracts array values between [ and ] (generally for yaml parsing)
///
/// \param Input string
/// \param delimiter
///
/// \return Vector of split double values
///
/////////////////////////////////////////////////////////////////////////////
std::vector<double> extract_list(std::string str);

/////////////////////////////////////////////////////////////////////////////
///
/// \fn ltrim
///
/// \brief <insert brief description>
///
/// <Insert longer more detailed description which
/// can span multiple lines if needed>
///
/// \param <function parameter description>
/// \param <function parameter description>
/// \param <function parameter description>
///
/// \return <return type and definition description if not void>
///
/////////////////////////////////////////////////////////////////////////////
std::string ltrim(const std::string& s);

/////////////////////////////////////////////////////////////////////////////
///
/// \fn rtrim
///
/// \brief <insert brief description>
///
/// <Insert longer more detailed description which
/// can span multiple lines if needed>
///
/// \param <function parameter description>
/// \param <function parameter description>
/// \param <function parameter description>
///
/// \return <return type and definition description if not void>
///
/////////////////////////////////////////////////////////////////////////////
std::string rtrim(const std::string& s);

/////////////////////////////////////////////////////////////////////////////
///
/// \fn trim
///
/// \brief <insert brief description>
///
/// <Insert longer more detailed description which
/// can span multiple lines if needed>
///
/// \param <function parameter description>
/// \param <function parameter description>
/// \param <function parameter description>
///
/// \return <return type and definition description if not void>
///
/////////////////////////////////////////////////////////////////////////////
std::string trim(const std::string& s);

#endif