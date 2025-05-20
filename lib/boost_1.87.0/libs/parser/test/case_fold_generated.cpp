// Copyright (c) 2024 T. Zachary Laine
//
// Distributed under the Boost Software License, Version 1.0. (See accompanying
// file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)

// Warning: This header is auto-generated (see misc/generate_case_fold_tests.py).

#include <boost/parser/parser.hpp>

#include <boost/core/lightweight_test.hpp>


int main()
{

// Initial code points from CaseFolding.txt
char32_t const cps[] = {
    0x0041, 0x0042, 0x0043, 0x0044, 0x0045, 0x0046, 0x0047, 0x0048, 
    0x0049, 0x004A, 0x004B, 0x004C, 0x004D, 0x004E, 0x004F, 0x0050, 
    0x0051, 0x0052, 0x0053, 0x0054, 0x0055, 0x0056, 0x0057, 0x0058, 
    0x0059, 0x005A, 0x00B5, 0x00C0, 0x00C1, 0x00C2, 0x00C3, 0x00C4, 
    0x00C5, 0x00C6, 0x00C7, 0x00C8, 0x00C9, 0x00CA, 0x00CB, 0x00CC, 
    0x00CD, 0x00CE, 0x00CF, 0x00D0, 0x00D1, 0x00D2, 0x00D3, 0x00D4, 
    0x00D5, 0x00D6, 0x00D8, 0x00D9, 0x00DA, 0x00DB, 0x00DC, 0x00DD, 
    0x00DE, 0x00DF, 0x0100, 0x0102, 0x0104, 0x0106, 0x0108, 0x010A, 
    0x010C, 0x010E, 0x0110, 0x0112, 0x0114, 0x0116, 0x0118, 0x011A, 
    0x011C, 0x011E, 0x0120, 0x0122, 0x0124, 0x0126, 0x0128, 0x012A, 
    0x012C, 0x012E, 0x0130, 0x0132, 0x0134, 0x0136, 0x0139, 0x013B, 
    0x013D, 0x013F, 0x0141, 0x0143, 0x0145, 0x0147, 0x0149, 0x014A, 
    0x014C, 0x014E, 0x0150, 0x0152, 0x0154, 0x0156, 0x0158, 0x015A, 
    0x015C, 0x015E, 0x0160, 0x0162, 0x0164, 0x0166, 0x0168, 0x016A, 
    0x016C, 0x016E, 0x0170, 0x0172, 0x0174, 0x0176, 0x0178, 0x0179, 
    0x017B, 0x017D, 0x017F, 0x0181, 0x0182, 0x0184, 0x0186, 0x0187, 
    0x0189, 0x018A, 0x018B, 0x018E, 0x018F, 0x0190, 0x0191, 0x0193, 
    0x0194, 0x0196, 0x0197, 0x0198, 0x019C, 0x019D, 0x019F, 0x01A0, 
    0x01A2, 0x01A4, 0x01A6, 0x01A7, 0x01A9, 0x01AC, 0x01AE, 0x01AF, 
    0x01B1, 0x01B2, 0x01B3, 0x01B5, 0x01B7, 0x01B8, 0x01BC, 0x01C4, 
    0x01C5, 0x01C7, 0x01C8, 0x01CA, 0x01CB, 0x01CD, 0x01CF, 0x01D1, 
    0x01D3, 0x01D5, 0x01D7, 0x01D9, 0x01DB, 0x01DE, 0x01E0, 0x01E2, 
    0x01E4, 0x01E6, 0x01E8, 0x01EA, 0x01EC, 0x01EE, 0x01F0, 0x01F1, 
    0x01F2, 0x01F4, 0x01F6, 0x01F7, 0x01F8, 0x01FA, 0x01FC, 0x01FE, 
    0x0200, 0x0202, 0x0204, 0x0206, 0x0208, 0x020A, 0x020C, 0x020E, 
    0x0210, 0x0212, 0x0214, 0x0216, 0x0218, 0x021A, 0x021C, 0x021E, 
    0x0220, 0x0222, 0x0224, 0x0226, 0x0228, 0x022A, 0x022C, 0x022E, 
    0x0230, 0x0232, 0x023A, 0x023B, 0x023D, 0x023E, 0x0241, 0x0243, 
    0x0244, 0x0245, 0x0246, 0x0248, 0x024A, 0x024C, 0x024E, 0x0345, 
    0x0370, 0x0372, 0x0376, 0x037F, 0x0386, 0x0388, 0x0389, 0x038A, 
    0x038C, 0x038E, 0x038F, 0x0390, 0x0391, 0x0392, 0x0393, 0x0394, 
    0x0395, 0x0396, 0x0397, 0x0398, 0x0399, 0x039A, 0x039B, 0x039C, 
    0x039D, 0x039E, 0x039F, 0x03A0, 0x03A1, 0x03A3, 0x03A4, 0x03A5, 
    0x03A6, 0x03A7, 0x03A8, 0x03A9, 0x03AA, 0x03AB, 0x03B0, 0x03C2, 
    0x03CF, 0x03D0, 0x03D1, 0x03D5, 0x03D6, 0x03D8, 0x03DA, 0x03DC, 
    0x03DE, 0x03E0, 0x03E2, 0x03E4, 0x03E6, 0x03E8, 0x03EA, 0x03EC, 
    0x03EE, 0x03F0, 0x03F1, 0x03F4, 0x03F5, 0x03F7, 0x03F9, 0x03FA, 
    0x03FD, 0x03FE, 0x03FF, 0x0400, 0x0401, 0x0402, 0x0403, 0x0404, 
    0x0405, 0x0406, 0x0407, 0x0408, 0x0409, 0x040A, 0x040B, 0x040C, 
    0x040D, 0x040E, 0x040F, 0x0410, 0x0411, 0x0412, 0x0413, 0x0414, 
    0x0415, 0x0416, 0x0417, 0x0418, 0x0419, 0x041A, 0x041B, 0x041C, 
    0x041D, 0x041E, 0x041F, 0x0420, 0x0421, 0x0422, 0x0423, 0x0424, 
    0x0425, 0x0426, 0x0427, 0x0428, 0x0429, 0x042A, 0x042B, 0x042C, 
    0x042D, 0x042E, 0x042F, 0x0460, 0x0462, 0x0464, 0x0466, 0x0468, 
    0x046A, 0x046C, 0x046E, 0x0470, 0x0472, 0x0474, 0x0476, 0x0478, 
    0x047A, 0x047C, 0x047E, 0x0480, 0x048A, 0x048C, 0x048E, 0x0490, 
    0x0492, 0x0494, 0x0496, 0x0498, 0x049A, 0x049C, 0x049E, 0x04A0, 
    0x04A2, 0x04A4, 0x04A6, 0x04A8, 0x04AA, 0x04AC, 0x04AE, 0x04B0, 
    0x04B2, 0x04B4, 0x04B6, 0x04B8, 0x04BA, 0x04BC, 0x04BE, 0x04C0, 
    0x04C1, 0x04C3, 0x04C5, 0x04C7, 0x04C9, 0x04CB, 0x04CD, 0x04D0, 
    0x04D2, 0x04D4, 0x04D6, 0x04D8, 0x04DA, 0x04DC, 0x04DE, 0x04E0, 
    0x04E2, 0x04E4, 0x04E6, 0x04E8, 0x04EA, 0x04EC, 0x04EE, 0x04F0, 
    0x04F2, 0x04F4, 0x04F6, 0x04F8, 0x04FA, 0x04FC, 0x04FE, 0x0500, 
    0x0502, 0x0504, 0x0506, 0x0508, 0x050A, 0x050C, 0x050E, 0x0510, 
    0x0512, 0x0514, 0x0516, 0x0518, 0x051A, 0x051C, 0x051E, 0x0520, 
    0x0522, 0x0524, 0x0526, 0x0528, 0x052A, 0x052C, 0x052E, 0x0531, 
    0x0532, 0x0533, 0x0534, 0x0535, 0x0536, 0x0537, 0x0538, 0x0539, 
    0x053A, 0x053B, 0x053C, 0x053D, 0x053E, 0x053F, 0x0540, 0x0541, 
    0x0542, 0x0543, 0x0544, 0x0545, 0x0546, 0x0547, 0x0548, 0x0549, 
    0x054A, 0x054B, 0x054C, 0x054D, 0x054E, 0x054F, 0x0550, 0x0551, 
    0x0552, 0x0553, 0x0554, 0x0555, 0x0556, 0x0587, 0x10A0, 0x10A1, 
    0x10A2, 0x10A3, 0x10A4, 0x10A5, 0x10A6, 0x10A7, 0x10A8, 0x10A9, 
    0x10AA, 0x10AB, 0x10AC, 0x10AD, 0x10AE, 0x10AF, 0x10B0, 0x10B1, 
    0x10B2, 0x10B3, 0x10B4, 0x10B5, 0x10B6, 0x10B7, 0x10B8, 0x10B9, 
    0x10BA, 0x10BB, 0x10BC, 0x10BD, 0x10BE, 0x10BF, 0x10C0, 0x10C1, 
    0x10C2, 0x10C3, 0x10C4, 0x10C5, 0x10C7, 0x10CD, 0x13F8, 0x13F9, 
    0x13FA, 0x13FB, 0x13FC, 0x13FD, 0x1C80, 0x1C81, 0x1C82, 0x1C83, 
    0x1C84, 0x1C85, 0x1C86, 0x1C87, 0x1C88, 0x1C90, 0x1C91, 0x1C92, 
    0x1C93, 0x1C94, 0x1C95, 0x1C96, 0x1C97, 0x1C98, 0x1C99, 0x1C9A, 
    0x1C9B, 0x1C9C, 0x1C9D, 0x1C9E, 0x1C9F, 0x1CA0, 0x1CA1, 0x1CA2, 
    0x1CA3, 0x1CA4, 0x1CA5, 0x1CA6, 0x1CA7, 0x1CA8, 0x1CA9, 0x1CAA, 
    0x1CAB, 0x1CAC, 0x1CAD, 0x1CAE, 0x1CAF, 0x1CB0, 0x1CB1, 0x1CB2, 
    0x1CB3, 0x1CB4, 0x1CB5, 0x1CB6, 0x1CB7, 0x1CB8, 0x1CB9, 0x1CBA, 
    0x1CBD, 0x1CBE, 0x1CBF, 0x1E00, 0x1E02, 0x1E04, 0x1E06, 0x1E08, 
    0x1E0A, 0x1E0C, 0x1E0E, 0x1E10, 0x1E12, 0x1E14, 0x1E16, 0x1E18, 
    0x1E1A, 0x1E1C, 0x1E1E, 0x1E20, 0x1E22, 0x1E24, 0x1E26, 0x1E28, 
    0x1E2A, 0x1E2C, 0x1E2E, 0x1E30, 0x1E32, 0x1E34, 0x1E36, 0x1E38, 
    0x1E3A, 0x1E3C, 0x1E3E, 0x1E40, 0x1E42, 0x1E44, 0x1E46, 0x1E48, 
    0x1E4A, 0x1E4C, 0x1E4E, 0x1E50, 0x1E52, 0x1E54, 0x1E56, 0x1E58, 
    0x1E5A, 0x1E5C, 0x1E5E, 0x1E60, 0x1E62, 0x1E64, 0x1E66, 0x1E68, 
    0x1E6A, 0x1E6C, 0x1E6E, 0x1E70, 0x1E72, 0x1E74, 0x1E76, 0x1E78, 
    0x1E7A, 0x1E7C, 0x1E7E, 0x1E80, 0x1E82, 0x1E84, 0x1E86, 0x1E88, 
    0x1E8A, 0x1E8C, 0x1E8E, 0x1E90, 0x1E92, 0x1E94, 0x1E96, 0x1E97, 
    0x1E98, 0x1E99, 0x1E9A, 0x1E9B, 0x1E9E, 0x1EA0, 0x1EA2, 0x1EA4, 
    0x1EA6, 0x1EA8, 0x1EAA, 0x1EAC, 0x1EAE, 0x1EB0, 0x1EB2, 0x1EB4, 
    0x1EB6, 0x1EB8, 0x1EBA, 0x1EBC, 0x1EBE, 0x1EC0, 0x1EC2, 0x1EC4, 
    0x1EC6, 0x1EC8, 0x1ECA, 0x1ECC, 0x1ECE, 0x1ED0, 0x1ED2, 0x1ED4, 
    0x1ED6, 0x1ED8, 0x1EDA, 0x1EDC, 0x1EDE, 0x1EE0, 0x1EE2, 0x1EE4, 
    0x1EE6, 0x1EE8, 0x1EEA, 0x1EEC, 0x1EEE, 0x1EF0, 0x1EF2, 0x1EF4, 
    0x1EF6, 0x1EF8, 0x1EFA, 0x1EFC, 0x1EFE, 0x1F08, 0x1F09, 0x1F0A, 
    0x1F0B, 0x1F0C, 0x1F0D, 0x1F0E, 0x1F0F, 0x1F18, 0x1F19, 0x1F1A, 
    0x1F1B, 0x1F1C, 0x1F1D, 0x1F28, 0x1F29, 0x1F2A, 0x1F2B, 0x1F2C, 
    0x1F2D, 0x1F2E, 0x1F2F, 0x1F38, 0x1F39, 0x1F3A, 0x1F3B, 0x1F3C, 
    0x1F3D, 0x1F3E, 0x1F3F, 0x1F48, 0x1F49, 0x1F4A, 0x1F4B, 0x1F4C, 
    0x1F4D, 0x1F50, 0x1F52, 0x1F54, 0x1F56, 0x1F59, 0x1F5B, 0x1F5D, 
    0x1F5F, 0x1F68, 0x1F69, 0x1F6A, 0x1F6B, 0x1F6C, 0x1F6D, 0x1F6E, 
    0x1F6F, 0x1F80, 0x1F81, 0x1F82, 0x1F83, 0x1F84, 0x1F85, 0x1F86, 
    0x1F87, 0x1F88, 0x1F89, 0x1F8A, 0x1F8B, 0x1F8C, 0x1F8D, 0x1F8E, 
    0x1F8F, 0x1F90, 0x1F91, 0x1F92, 0x1F93, 0x1F94, 0x1F95, 0x1F96, 
    0x1F97, 0x1F98, 0x1F99, 0x1F9A, 0x1F9B, 0x1F9C, 0x1F9D, 0x1F9E, 
    0x1F9F, 0x1FA0, 0x1FA1, 0x1FA2, 0x1FA3, 0x1FA4, 0x1FA5, 0x1FA6, 
    0x1FA7, 0x1FA8, 0x1FA9, 0x1FAA, 0x1FAB, 0x1FAC, 0x1FAD, 0x1FAE, 
    0x1FAF, 0x1FB2, 0x1FB3, 0x1FB4, 0x1FB6, 0x1FB7, 0x1FB8, 0x1FB9, 
    0x1FBA, 0x1FBB, 0x1FBC, 0x1FBE, 0x1FC2, 0x1FC3, 0x1FC4, 0x1FC6, 
    0x1FC7, 0x1FC8, 0x1FC9, 0x1FCA, 0x1FCB, 0x1FCC, 0x1FD2, 0x1FD3, 
    0x1FD6, 0x1FD7, 0x1FD8, 0x1FD9, 0x1FDA, 0x1FDB, 0x1FE2, 0x1FE3, 
    0x1FE4, 0x1FE6, 0x1FE7, 0x1FE8, 0x1FE9, 0x1FEA, 0x1FEB, 0x1FEC, 
    0x1FF2, 0x1FF3, 0x1FF4, 0x1FF6, 0x1FF7, 0x1FF8, 0x1FF9, 0x1FFA, 
    0x1FFB, 0x1FFC, 0x2126, 0x212A, 0x212B, 0x2132, 0x2160, 0x2161, 
    0x2162, 0x2163, 0x2164, 0x2165, 0x2166, 0x2167, 0x2168, 0x2169, 
    0x216A, 0x216B, 0x216C, 0x216D, 0x216E, 0x216F, 0x2183, 0x24B6, 
    0x24B7, 0x24B8, 0x24B9, 0x24BA, 0x24BB, 0x24BC, 0x24BD, 0x24BE, 
    0x24BF, 0x24C0, 0x24C1, 0x24C2, 0x24C3, 0x24C4, 0x24C5, 0x24C6, 
    0x24C7, 0x24C8, 0x24C9, 0x24CA, 0x24CB, 0x24CC, 0x24CD, 0x24CE, 
    0x24CF, 0x2C00, 0x2C01, 0x2C02, 0x2C03, 0x2C04, 0x2C05, 0x2C06, 
    0x2C07, 0x2C08, 0x2C09, 0x2C0A, 0x2C0B, 0x2C0C, 0x2C0D, 0x2C0E, 
    0x2C0F, 0x2C10, 0x2C11, 0x2C12, 0x2C13, 0x2C14, 0x2C15, 0x2C16, 
    0x2C17, 0x2C18, 0x2C19, 0x2C1A, 0x2C1B, 0x2C1C, 0x2C1D, 0x2C1E, 
    0x2C1F, 0x2C20, 0x2C21, 0x2C22, 0x2C23, 0x2C24, 0x2C25, 0x2C26, 
    0x2C27, 0x2C28, 0x2C29, 0x2C2A, 0x2C2B, 0x2C2C, 0x2C2D, 0x2C2E, 
    0x2C2F, 0x2C60, 0x2C62, 0x2C63, 0x2C64, 0x2C67, 0x2C69, 0x2C6B, 
    0x2C6D, 0x2C6E, 0x2C6F, 0x2C70, 0x2C72, 0x2C75, 0x2C7E, 0x2C7F, 
    0x2C80, 0x2C82, 0x2C84, 0x2C86, 0x2C88, 0x2C8A, 0x2C8C, 0x2C8E, 
    0x2C90, 0x2C92, 0x2C94, 0x2C96, 0x2C98, 0x2C9A, 0x2C9C, 0x2C9E, 
    0x2CA0, 0x2CA2, 0x2CA4, 0x2CA6, 0x2CA8, 0x2CAA, 0x2CAC, 0x2CAE, 
    0x2CB0, 0x2CB2, 0x2CB4, 0x2CB6, 0x2CB8, 0x2CBA, 0x2CBC, 0x2CBE, 
    0x2CC0, 0x2CC2, 0x2CC4, 0x2CC6, 0x2CC8, 0x2CCA, 0x2CCC, 0x2CCE, 
    0x2CD0, 0x2CD2, 0x2CD4, 0x2CD6, 0x2CD8, 0x2CDA, 0x2CDC, 0x2CDE, 
    0x2CE0, 0x2CE2, 0x2CEB, 0x2CED, 0x2CF2, 0xA640, 0xA642, 0xA644, 
    0xA646, 0xA648, 0xA64A, 0xA64C, 0xA64E, 0xA650, 0xA652, 0xA654, 
    0xA656, 0xA658, 0xA65A, 0xA65C, 0xA65E, 0xA660, 0xA662, 0xA664, 
    0xA666, 0xA668, 0xA66A, 0xA66C, 0xA680, 0xA682, 0xA684, 0xA686, 
    0xA688, 0xA68A, 0xA68C, 0xA68E, 0xA690, 0xA692, 0xA694, 0xA696, 
    0xA698, 0xA69A, 0xA722, 0xA724, 0xA726, 0xA728, 0xA72A, 0xA72C, 
    0xA72E, 0xA732, 0xA734, 0xA736, 0xA738, 0xA73A, 0xA73C, 0xA73E, 
    0xA740, 0xA742, 0xA744, 0xA746, 0xA748, 0xA74A, 0xA74C, 0xA74E, 
    0xA750, 0xA752, 0xA754, 0xA756, 0xA758, 0xA75A, 0xA75C, 0xA75E, 
    0xA760, 0xA762, 0xA764, 0xA766, 0xA768, 0xA76A, 0xA76C, 0xA76E, 
    0xA779, 0xA77B, 0xA77D, 0xA77E, 0xA780, 0xA782, 0xA784, 0xA786, 
    0xA78B, 0xA78D, 0xA790, 0xA792, 0xA796, 0xA798, 0xA79A, 0xA79C, 
    0xA79E, 0xA7A0, 0xA7A2, 0xA7A4, 0xA7A6, 0xA7A8, 0xA7AA, 0xA7AB, 
    0xA7AC, 0xA7AD, 0xA7AE, 0xA7B0, 0xA7B1, 0xA7B2, 0xA7B3, 0xA7B4, 
    0xA7B6, 0xA7B8, 0xA7BA, 0xA7BC, 0xA7BE, 0xA7C0, 0xA7C2, 0xA7C4, 
    0xA7C5, 0xA7C6, 0xA7C7, 0xA7C9, 0xA7D0, 0xA7D6, 0xA7D8, 0xA7F5, 
    0xAB70, 0xAB71, 0xAB72, 0xAB73, 0xAB74, 0xAB75, 0xAB76, 0xAB77, 
    0xAB78, 0xAB79, 0xAB7A, 0xAB7B, 0xAB7C, 0xAB7D, 0xAB7E, 0xAB7F, 
    0xAB80, 0xAB81, 0xAB82, 0xAB83, 0xAB84, 0xAB85, 0xAB86, 0xAB87, 
    0xAB88, 0xAB89, 0xAB8A, 0xAB8B, 0xAB8C, 0xAB8D, 0xAB8E, 0xAB8F, 
    0xAB90, 0xAB91, 0xAB92, 0xAB93, 0xAB94, 0xAB95, 0xAB96, 0xAB97, 
    0xAB98, 0xAB99, 0xAB9A, 0xAB9B, 0xAB9C, 0xAB9D, 0xAB9E, 0xAB9F, 
    0xABA0, 0xABA1, 0xABA2, 0xABA3, 0xABA4, 0xABA5, 0xABA6, 0xABA7, 
    0xABA8, 0xABA9, 0xABAA, 0xABAB, 0xABAC, 0xABAD, 0xABAE, 0xABAF, 
    0xABB0, 0xABB1, 0xABB2, 0xABB3, 0xABB4, 0xABB5, 0xABB6, 0xABB7, 
    0xABB8, 0xABB9, 0xABBA, 0xABBB, 0xABBC, 0xABBD, 0xABBE, 0xABBF, 
    0xFB00, 0xFB01, 0xFB02, 0xFB03, 0xFB04, 0xFB05, 0xFB06, 0xFB13, 
    0xFB14, 0xFB15, 0xFB16, 0xFB17, 0xFF21, 0xFF22, 0xFF23, 0xFF24, 
    0xFF25, 0xFF26, 0xFF27, 0xFF28, 0xFF29, 0xFF2A, 0xFF2B, 0xFF2C, 
    0xFF2D, 0xFF2E, 0xFF2F, 0xFF30, 0xFF31, 0xFF32, 0xFF33, 0xFF34, 
    0xFF35, 0xFF36, 0xFF37, 0xFF38, 0xFF39, 0xFF3A, 0x10400, 0x10401, 
    0x10402, 0x10403, 0x10404, 0x10405, 0x10406, 0x10407, 0x10408, 0x10409, 
    0x1040A, 0x1040B, 0x1040C, 0x1040D, 0x1040E, 0x1040F, 0x10410, 0x10411, 
    0x10412, 0x10413, 0x10414, 0x10415, 0x10416, 0x10417, 0x10418, 0x10419, 
    0x1041A, 0x1041B, 0x1041C, 0x1041D, 0x1041E, 0x1041F, 0x10420, 0x10421, 
    0x10422, 0x10423, 0x10424, 0x10425, 0x10426, 0x10427, 0x104B0, 0x104B1, 
    0x104B2, 0x104B3, 0x104B4, 0x104B5, 0x104B6, 0x104B7, 0x104B8, 0x104B9, 
    0x104BA, 0x104BB, 0x104BC, 0x104BD, 0x104BE, 0x104BF, 0x104C0, 0x104C1, 
    0x104C2, 0x104C3, 0x104C4, 0x104C5, 0x104C6, 0x104C7, 0x104C8, 0x104C9, 
    0x104CA, 0x104CB, 0x104CC, 0x104CD, 0x104CE, 0x104CF, 0x104D0, 0x104D1, 
    0x104D2, 0x104D3, 0x10570, 0x10571, 0x10572, 0x10573, 0x10574, 0x10575, 
    0x10576, 0x10577, 0x10578, 0x10579, 0x1057A, 0x1057C, 0x1057D, 0x1057E, 
    0x1057F, 0x10580, 0x10581, 0x10582, 0x10583, 0x10584, 0x10585, 0x10586, 
    0x10587, 0x10588, 0x10589, 0x1058A, 0x1058C, 0x1058D, 0x1058E, 0x1058F, 
    0x10590, 0x10591, 0x10592, 0x10594, 0x10595, 0x10C80, 0x10C81, 0x10C82, 
    0x10C83, 0x10C84, 0x10C85, 0x10C86, 0x10C87, 0x10C88, 0x10C89, 0x10C8A, 
    0x10C8B, 0x10C8C, 0x10C8D, 0x10C8E, 0x10C8F, 0x10C90, 0x10C91, 0x10C92, 
    0x10C93, 0x10C94, 0x10C95, 0x10C96, 0x10C97, 0x10C98, 0x10C99, 0x10C9A, 
    0x10C9B, 0x10C9C, 0x10C9D, 0x10C9E, 0x10C9F, 0x10CA0, 0x10CA1, 0x10CA2, 
    0x10CA3, 0x10CA4, 0x10CA5, 0x10CA6, 0x10CA7, 0x10CA8, 0x10CA9, 0x10CAA, 
    0x10CAB, 0x10CAC, 0x10CAD, 0x10CAE, 0x10CAF, 0x10CB0, 0x10CB1, 0x10CB2, 
    0x118A0, 0x118A1, 0x118A2, 0x118A3, 0x118A4, 0x118A5, 0x118A6, 0x118A7, 
    0x118A8, 0x118A9, 0x118AA, 0x118AB, 0x118AC, 0x118AD, 0x118AE, 0x118AF, 
    0x118B0, 0x118B1, 0x118B2, 0x118B3, 0x118B4, 0x118B5, 0x118B6, 0x118B7, 
    0x118B8, 0x118B9, 0x118BA, 0x118BB, 0x118BC, 0x118BD, 0x118BE, 0x118BF, 
    0x16E40, 0x16E41, 0x16E42, 0x16E43, 0x16E44, 0x16E45, 0x16E46, 0x16E47, 
    0x16E48, 0x16E49, 0x16E4A, 0x16E4B, 0x16E4C, 0x16E4D, 0x16E4E, 0x16E4F, 
    0x16E50, 0x16E51, 0x16E52, 0x16E53, 0x16E54, 0x16E55, 0x16E56, 0x16E57, 
    0x16E58, 0x16E59, 0x16E5A, 0x16E5B, 0x16E5C, 0x16E5D, 0x16E5E, 0x16E5F, 
    0x1E900, 0x1E901, 0x1E902, 0x1E903, 0x1E904, 0x1E905, 0x1E906, 0x1E907, 
    0x1E908, 0x1E909, 0x1E90A, 0x1E90B, 0x1E90C, 0x1E90D, 0x1E90E, 0x1E90F, 
    0x1E910, 0x1E911, 0x1E912, 0x1E913, 0x1E914, 0x1E915, 0x1E916, 0x1E917, 
    0x1E918, 0x1E919, 0x1E91A, 0x1E91B, 0x1E91C, 0x1E91D, 0x1E91E, 0x1E91F, 
    0x1E920, 0x1E921, 
};

[[maybe_unused]] char32_t const max_test_cp = 0x1E921 + 100;

// hits_0)
{
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x0061, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x0041, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x0062, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x0042, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x0063, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x0043, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x0064, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x0044, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x0065, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x0045, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x0066, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x0046, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x0067, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x0047, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x0068, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x0048, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x0069, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x0049, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x006A, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x004A, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x006B, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x004B, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x006C, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x004C, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x006D, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x004D, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x006E, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x004E, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x006F, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x004F, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x0070, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x0050, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x0071, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x0051, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x0072, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x0052, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x0073, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x0053, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x0074, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x0054, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x0075, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x0055, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x0076, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x0056, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x0077, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x0057, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x0078, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x0058, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x0079, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x0059, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x007A, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x005A, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x03BC, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x00B5, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x00E0, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x00C0, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x00E1, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x00C1, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x00E2, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x00C2, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x00E3, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x00C3, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x00E4, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x00C4, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x00E5, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x00C5, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x00E6, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x00C6, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x00E7, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x00C7, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x00E8, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x00C8, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x00E9, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x00C9, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x00EA, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x00CA, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x00EB, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x00CB, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x00EC, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x00CC, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x00ED, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x00CD, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x00EE, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x00CE, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x00EF, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x00CF, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x00F0, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x00D0, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x00F1, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x00D1, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x00F2, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x00D2, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x00F3, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x00D3, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x00F4, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x00D4, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x00F5, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x00D5, result.begin());
        BOOST_TEST(result == expected);
    }
}

// test_1)
{
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x00F6, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x00D6, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x00F8, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x00D8, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x00F9, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x00D9, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x00FA, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x00DA, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x00FB, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x00DB, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x00FC, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x00DC, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x00FD, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x00DD, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x00FE, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x00DE, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x0073, 0x0073, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x00DF, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x0101, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x0100, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x0103, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x0102, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x0105, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x0104, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x0107, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x0106, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x0109, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x0108, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x010B, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x010A, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x010D, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x010C, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x010F, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x010E, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x0111, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x0110, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x0113, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x0112, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x0115, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x0114, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x0117, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x0116, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x0119, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x0118, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x011B, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x011A, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x011D, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x011C, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x011F, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x011E, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x0121, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x0120, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x0123, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x0122, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x0125, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x0124, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x0127, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x0126, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x0129, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x0128, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x012B, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x012A, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x012D, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x012C, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x012F, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x012E, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x0069, 0x0307, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x0130, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x0133, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x0132, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x0135, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x0134, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x0137, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x0136, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x013A, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x0139, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x013C, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x013B, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x013E, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x013D, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x0140, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x013F, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x0142, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x0141, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x0144, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x0143, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x0146, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x0145, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x0148, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x0147, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x02BC, 0x006E, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x0149, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x014B, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x014A, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x014D, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x014C, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x014F, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x014E, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x0151, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x0150, result.begin());
        BOOST_TEST(result == expected);
    }
}

// test_2)
{
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x0153, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x0152, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x0155, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x0154, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x0157, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x0156, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x0159, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x0158, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x015B, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x015A, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x015D, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x015C, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x015F, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x015E, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x0161, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x0160, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x0163, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x0162, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x0165, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x0164, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x0167, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x0166, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x0169, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x0168, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x016B, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x016A, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x016D, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x016C, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x016F, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x016E, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x0171, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x0170, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x0173, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x0172, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x0175, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x0174, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x0177, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x0176, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x00FF, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x0178, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x017A, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x0179, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x017C, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x017B, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x017E, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x017D, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x0073, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x017F, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x0253, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x0181, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x0183, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x0182, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x0185, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x0184, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x0254, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x0186, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x0188, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x0187, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x0256, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x0189, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x0257, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x018A, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x018C, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x018B, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x01DD, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x018E, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x0259, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x018F, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x025B, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x0190, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x0192, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x0191, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x0260, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x0193, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x0263, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x0194, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x0269, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x0196, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x0268, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x0197, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x0199, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x0198, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x026F, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x019C, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x0272, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x019D, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x0275, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x019F, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x01A1, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x01A0, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x01A3, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x01A2, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x01A5, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x01A4, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x0280, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x01A6, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x01A8, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x01A7, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x0283, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x01A9, result.begin());
        BOOST_TEST(result == expected);
    }
}

// test_3)
{
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x01AD, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x01AC, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x0288, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x01AE, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x01B0, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x01AF, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x028A, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x01B1, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x028B, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x01B2, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x01B4, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x01B3, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x01B6, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x01B5, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x0292, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x01B7, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x01B9, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x01B8, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x01BD, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x01BC, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x01C6, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x01C4, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x01C6, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x01C5, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x01C9, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x01C7, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x01C9, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x01C8, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x01CC, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x01CA, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x01CC, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x01CB, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x01CE, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x01CD, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x01D0, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x01CF, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x01D2, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x01D1, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x01D4, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x01D3, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x01D6, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x01D5, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x01D8, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x01D7, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x01DA, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x01D9, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x01DC, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x01DB, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x01DF, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x01DE, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x01E1, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x01E0, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x01E3, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x01E2, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x01E5, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x01E4, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x01E7, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x01E6, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x01E9, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x01E8, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x01EB, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x01EA, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x01ED, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x01EC, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x01EF, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x01EE, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x006A, 0x030C, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x01F0, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x01F3, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x01F1, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x01F3, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x01F2, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x01F5, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x01F4, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x0195, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x01F6, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x01BF, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x01F7, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x01F9, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x01F8, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x01FB, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x01FA, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x01FD, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x01FC, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x01FF, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x01FE, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x0201, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x0200, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x0203, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x0202, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x0205, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x0204, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x0207, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x0206, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x0209, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x0208, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x020B, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x020A, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x020D, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x020C, result.begin());
        BOOST_TEST(result == expected);
    }
}

// test_4)
{
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x020F, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x020E, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x0211, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x0210, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x0213, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x0212, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x0215, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x0214, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x0217, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x0216, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x0219, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x0218, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x021B, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x021A, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x021D, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x021C, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x021F, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x021E, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x019E, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x0220, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x0223, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x0222, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x0225, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x0224, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x0227, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x0226, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x0229, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x0228, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x022B, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x022A, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x022D, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x022C, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x022F, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x022E, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x0231, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x0230, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x0233, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x0232, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x2C65, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x023A, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x023C, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x023B, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x019A, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x023D, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x2C66, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x023E, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x0242, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x0241, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x0180, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x0243, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x0289, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x0244, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x028C, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x0245, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x0247, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x0246, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x0249, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x0248, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x024B, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x024A, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x024D, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x024C, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x024F, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x024E, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x03B9, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x0345, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x0371, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x0370, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x0373, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x0372, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x0377, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x0376, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x03F3, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x037F, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x03AC, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x0386, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x03AD, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x0388, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x03AE, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x0389, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x03AF, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x038A, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x03CC, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x038C, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x03CD, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x038E, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x03CE, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x038F, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x03B9, 0x0308, 0x0301, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x0390, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x03B1, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x0391, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x03B2, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x0392, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x03B3, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x0393, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x03B4, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x0394, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x03B5, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x0395, result.begin());
        BOOST_TEST(result == expected);
    }
}

// test_5)
{
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x03B6, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x0396, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x03B7, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x0397, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x03B8, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x0398, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x03B9, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x0399, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x03BA, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x039A, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x03BB, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x039B, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x03BC, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x039C, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x03BD, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x039D, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x03BE, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x039E, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x03BF, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x039F, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x03C0, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x03A0, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x03C1, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x03A1, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x03C3, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x03A3, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x03C4, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x03A4, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x03C5, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x03A5, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x03C6, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x03A6, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x03C7, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x03A7, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x03C8, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x03A8, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x03C9, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x03A9, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x03CA, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x03AA, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x03CB, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x03AB, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x03C5, 0x0308, 0x0301, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x03B0, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x03C3, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x03C2, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x03D7, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x03CF, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x03B2, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x03D0, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x03B8, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x03D1, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x03C6, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x03D5, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x03C0, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x03D6, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x03D9, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x03D8, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x03DB, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x03DA, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x03DD, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x03DC, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x03DF, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x03DE, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x03E1, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x03E0, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x03E3, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x03E2, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x03E5, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x03E4, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x03E7, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x03E6, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x03E9, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x03E8, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x03EB, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x03EA, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x03ED, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x03EC, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x03EF, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x03EE, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x03BA, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x03F0, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x03C1, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x03F1, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x03B8, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x03F4, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x03B5, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x03F5, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x03F8, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x03F7, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x03F2, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x03F9, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x03FB, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x03FA, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x037B, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x03FD, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x037C, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x03FE, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x037D, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x03FF, result.begin());
        BOOST_TEST(result == expected);
    }
}

// test_6)
{
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x0450, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x0400, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x0451, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x0401, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x0452, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x0402, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x0453, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x0403, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x0454, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x0404, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x0455, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x0405, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x0456, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x0406, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x0457, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x0407, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x0458, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x0408, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x0459, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x0409, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x045A, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x040A, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x045B, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x040B, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x045C, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x040C, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x045D, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x040D, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x045E, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x040E, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x045F, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x040F, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x0430, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x0410, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x0431, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x0411, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x0432, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x0412, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x0433, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x0413, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x0434, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x0414, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x0435, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x0415, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x0436, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x0416, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x0437, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x0417, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x0438, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x0418, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x0439, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x0419, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x043A, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x041A, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x043B, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x041B, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x043C, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x041C, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x043D, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x041D, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x043E, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x041E, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x043F, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x041F, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x0440, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x0420, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x0441, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x0421, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x0442, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x0422, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x0443, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x0423, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x0444, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x0424, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x0445, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x0425, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x0446, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x0426, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x0447, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x0427, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x0448, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x0428, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x0449, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x0429, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x044A, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x042A, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x044B, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x042B, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x044C, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x042C, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x044D, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x042D, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x044E, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x042E, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x044F, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x042F, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x0461, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x0460, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x0463, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x0462, result.begin());
        BOOST_TEST(result == expected);
    }
}

// test_7)
{
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x0465, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x0464, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x0467, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x0466, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x0469, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x0468, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x046B, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x046A, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x046D, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x046C, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x046F, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x046E, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x0471, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x0470, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x0473, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x0472, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x0475, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x0474, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x0477, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x0476, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x0479, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x0478, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x047B, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x047A, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x047D, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x047C, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x047F, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x047E, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x0481, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x0480, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x048B, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x048A, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x048D, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x048C, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x048F, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x048E, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x0491, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x0490, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x0493, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x0492, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x0495, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x0494, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x0497, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x0496, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x0499, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x0498, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x049B, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x049A, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x049D, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x049C, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x049F, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x049E, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x04A1, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x04A0, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x04A3, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x04A2, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x04A5, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x04A4, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x04A7, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x04A6, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x04A9, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x04A8, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x04AB, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x04AA, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x04AD, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x04AC, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x04AF, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x04AE, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x04B1, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x04B0, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x04B3, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x04B2, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x04B5, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x04B4, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x04B7, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x04B6, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x04B9, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x04B8, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x04BB, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x04BA, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x04BD, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x04BC, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x04BF, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x04BE, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x04CF, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x04C0, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x04C2, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x04C1, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x04C4, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x04C3, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x04C6, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x04C5, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x04C8, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x04C7, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x04CA, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x04C9, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x04CC, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x04CB, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x04CE, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x04CD, result.begin());
        BOOST_TEST(result == expected);
    }
}

// test_8)
{
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x04D1, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x04D0, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x04D3, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x04D2, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x04D5, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x04D4, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x04D7, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x04D6, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x04D9, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x04D8, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x04DB, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x04DA, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x04DD, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x04DC, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x04DF, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x04DE, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x04E1, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x04E0, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x04E3, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x04E2, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x04E5, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x04E4, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x04E7, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x04E6, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x04E9, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x04E8, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x04EB, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x04EA, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x04ED, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x04EC, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x04EF, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x04EE, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x04F1, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x04F0, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x04F3, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x04F2, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x04F5, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x04F4, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x04F7, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x04F6, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x04F9, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x04F8, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x04FB, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x04FA, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x04FD, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x04FC, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x04FF, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x04FE, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x0501, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x0500, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x0503, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x0502, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x0505, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x0504, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x0507, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x0506, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x0509, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x0508, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x050B, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x050A, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x050D, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x050C, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x050F, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x050E, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x0511, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x0510, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x0513, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x0512, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x0515, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x0514, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x0517, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x0516, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x0519, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x0518, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x051B, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x051A, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x051D, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x051C, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x051F, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x051E, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x0521, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x0520, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x0523, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x0522, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x0525, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x0524, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x0527, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x0526, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x0529, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x0528, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x052B, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x052A, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x052D, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x052C, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x052F, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x052E, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x0561, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x0531, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x0562, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x0532, result.begin());
        BOOST_TEST(result == expected);
    }
}

// test_9)
{
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x0563, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x0533, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x0564, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x0534, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x0565, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x0535, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x0566, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x0536, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x0567, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x0537, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x0568, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x0538, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x0569, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x0539, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x056A, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x053A, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x056B, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x053B, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x056C, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x053C, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x056D, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x053D, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x056E, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x053E, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x056F, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x053F, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x0570, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x0540, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x0571, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x0541, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x0572, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x0542, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x0573, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x0543, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x0574, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x0544, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x0575, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x0545, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x0576, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x0546, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x0577, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x0547, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x0578, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x0548, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x0579, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x0549, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x057A, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x054A, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x057B, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x054B, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x057C, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x054C, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x057D, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x054D, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x057E, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x054E, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x057F, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x054F, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x0580, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x0550, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x0581, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x0551, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x0582, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x0552, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x0583, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x0553, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x0584, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x0554, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x0585, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x0555, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x0586, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x0556, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x0565, 0x0582, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x0587, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x2D00, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x10A0, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x2D01, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x10A1, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x2D02, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x10A2, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x2D03, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x10A3, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x2D04, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x10A4, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x2D05, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x10A5, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x2D06, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x10A6, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x2D07, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x10A7, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x2D08, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x10A8, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x2D09, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x10A9, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x2D0A, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x10AA, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x2D0B, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x10AB, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x2D0C, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x10AC, result.begin());
        BOOST_TEST(result == expected);
    }
}

// test_10)
{
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x2D0D, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x10AD, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x2D0E, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x10AE, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x2D0F, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x10AF, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x2D10, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x10B0, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x2D11, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x10B1, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x2D12, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x10B2, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x2D13, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x10B3, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x2D14, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x10B4, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x2D15, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x10B5, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x2D16, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x10B6, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x2D17, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x10B7, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x2D18, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x10B8, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x2D19, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x10B9, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x2D1A, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x10BA, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x2D1B, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x10BB, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x2D1C, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x10BC, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x2D1D, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x10BD, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x2D1E, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x10BE, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x2D1F, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x10BF, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x2D20, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x10C0, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x2D21, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x10C1, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x2D22, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x10C2, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x2D23, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x10C3, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x2D24, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x10C4, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x2D25, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x10C5, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x2D27, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x10C7, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x2D2D, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x10CD, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x13F0, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x13F8, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x13F1, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x13F9, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x13F2, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x13FA, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x13F3, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x13FB, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x13F4, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x13FC, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x13F5, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x13FD, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x0432, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1C80, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x0434, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1C81, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x043E, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1C82, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x0441, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1C83, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x0442, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1C84, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x0442, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1C85, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x044A, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1C86, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x0463, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1C87, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0xA64B, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1C88, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x10D0, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1C90, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x10D1, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1C91, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x10D2, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1C92, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x10D3, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1C93, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x10D4, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1C94, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x10D5, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1C95, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x10D6, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1C96, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x10D7, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1C97, result.begin());
        BOOST_TEST(result == expected);
    }
}

// test_11)
{
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x10D8, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1C98, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x10D9, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1C99, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x10DA, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1C9A, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x10DB, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1C9B, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x10DC, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1C9C, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x10DD, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1C9D, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x10DE, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1C9E, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x10DF, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1C9F, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x10E0, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1CA0, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x10E1, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1CA1, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x10E2, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1CA2, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x10E3, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1CA3, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x10E4, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1CA4, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x10E5, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1CA5, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x10E6, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1CA6, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x10E7, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1CA7, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x10E8, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1CA8, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x10E9, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1CA9, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x10EA, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1CAA, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x10EB, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1CAB, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x10EC, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1CAC, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x10ED, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1CAD, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x10EE, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1CAE, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x10EF, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1CAF, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x10F0, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1CB0, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x10F1, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1CB1, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x10F2, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1CB2, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x10F3, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1CB3, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x10F4, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1CB4, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x10F5, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1CB5, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x10F6, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1CB6, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x10F7, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1CB7, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x10F8, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1CB8, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x10F9, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1CB9, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x10FA, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1CBA, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x10FD, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1CBD, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x10FE, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1CBE, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x10FF, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1CBF, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1E01, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1E00, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1E03, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1E02, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1E05, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1E04, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1E07, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1E06, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1E09, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1E08, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1E0B, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1E0A, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1E0D, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1E0C, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1E0F, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1E0E, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1E11, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1E10, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1E13, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1E12, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1E15, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1E14, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1E17, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1E16, result.begin());
        BOOST_TEST(result == expected);
    }
}

// test_12)
{
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1E19, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1E18, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1E1B, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1E1A, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1E1D, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1E1C, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1E1F, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1E1E, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1E21, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1E20, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1E23, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1E22, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1E25, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1E24, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1E27, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1E26, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1E29, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1E28, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1E2B, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1E2A, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1E2D, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1E2C, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1E2F, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1E2E, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1E31, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1E30, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1E33, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1E32, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1E35, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1E34, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1E37, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1E36, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1E39, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1E38, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1E3B, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1E3A, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1E3D, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1E3C, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1E3F, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1E3E, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1E41, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1E40, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1E43, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1E42, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1E45, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1E44, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1E47, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1E46, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1E49, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1E48, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1E4B, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1E4A, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1E4D, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1E4C, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1E4F, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1E4E, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1E51, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1E50, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1E53, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1E52, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1E55, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1E54, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1E57, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1E56, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1E59, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1E58, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1E5B, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1E5A, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1E5D, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1E5C, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1E5F, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1E5E, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1E61, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1E60, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1E63, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1E62, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1E65, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1E64, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1E67, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1E66, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1E69, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1E68, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1E6B, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1E6A, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1E6D, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1E6C, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1E6F, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1E6E, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1E71, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1E70, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1E73, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1E72, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1E75, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1E74, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1E77, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1E76, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1E79, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1E78, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1E7B, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1E7A, result.begin());
        BOOST_TEST(result == expected);
    }
}

// test_13)
{
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1E7D, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1E7C, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1E7F, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1E7E, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1E81, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1E80, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1E83, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1E82, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1E85, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1E84, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1E87, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1E86, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1E89, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1E88, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1E8B, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1E8A, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1E8D, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1E8C, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1E8F, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1E8E, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1E91, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1E90, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1E93, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1E92, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1E95, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1E94, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x0068, 0x0331, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1E96, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x0074, 0x0308, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1E97, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x0077, 0x030A, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1E98, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x0079, 0x030A, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1E99, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x0061, 0x02BE, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1E9A, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1E61, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1E9B, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x0073, 0x0073, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1E9E, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1EA1, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1EA0, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1EA3, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1EA2, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1EA5, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1EA4, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1EA7, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1EA6, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1EA9, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1EA8, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1EAB, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1EAA, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1EAD, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1EAC, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1EAF, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1EAE, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1EB1, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1EB0, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1EB3, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1EB2, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1EB5, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1EB4, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1EB7, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1EB6, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1EB9, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1EB8, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1EBB, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1EBA, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1EBD, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1EBC, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1EBF, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1EBE, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1EC1, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1EC0, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1EC3, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1EC2, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1EC5, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1EC4, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1EC7, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1EC6, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1EC9, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1EC8, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1ECB, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1ECA, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1ECD, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1ECC, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1ECF, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1ECE, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1ED1, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1ED0, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1ED3, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1ED2, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1ED5, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1ED4, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1ED7, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1ED6, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1ED9, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1ED8, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1EDB, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1EDA, result.begin());
        BOOST_TEST(result == expected);
    }
}

// test_14)
{
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1EDD, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1EDC, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1EDF, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1EDE, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1EE1, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1EE0, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1EE3, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1EE2, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1EE5, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1EE4, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1EE7, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1EE6, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1EE9, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1EE8, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1EEB, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1EEA, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1EED, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1EEC, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1EEF, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1EEE, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1EF1, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1EF0, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1EF3, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1EF2, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1EF5, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1EF4, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1EF7, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1EF6, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1EF9, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1EF8, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1EFB, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1EFA, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1EFD, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1EFC, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1EFF, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1EFE, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1F00, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1F08, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1F01, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1F09, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1F02, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1F0A, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1F03, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1F0B, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1F04, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1F0C, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1F05, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1F0D, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1F06, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1F0E, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1F07, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1F0F, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1F10, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1F18, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1F11, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1F19, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1F12, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1F1A, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1F13, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1F1B, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1F14, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1F1C, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1F15, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1F1D, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1F20, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1F28, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1F21, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1F29, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1F22, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1F2A, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1F23, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1F2B, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1F24, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1F2C, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1F25, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1F2D, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1F26, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1F2E, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1F27, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1F2F, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1F30, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1F38, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1F31, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1F39, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1F32, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1F3A, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1F33, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1F3B, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1F34, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1F3C, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1F35, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1F3D, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1F36, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1F3E, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1F37, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1F3F, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1F40, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1F48, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1F41, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1F49, result.begin());
        BOOST_TEST(result == expected);
    }
}

// test_15)
{
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1F42, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1F4A, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1F43, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1F4B, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1F44, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1F4C, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1F45, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1F4D, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x03C5, 0x0313, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1F50, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x03C5, 0x0313, 0x0300, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1F52, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x03C5, 0x0313, 0x0301, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1F54, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x03C5, 0x0313, 0x0342, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1F56, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1F51, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1F59, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1F53, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1F5B, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1F55, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1F5D, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1F57, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1F5F, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1F60, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1F68, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1F61, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1F69, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1F62, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1F6A, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1F63, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1F6B, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1F64, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1F6C, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1F65, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1F6D, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1F66, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1F6E, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1F67, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1F6F, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1F00, 0x03B9, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1F80, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1F01, 0x03B9, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1F81, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1F02, 0x03B9, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1F82, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1F03, 0x03B9, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1F83, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1F04, 0x03B9, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1F84, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1F05, 0x03B9, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1F85, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1F06, 0x03B9, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1F86, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1F07, 0x03B9, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1F87, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1F00, 0x03B9, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1F88, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1F01, 0x03B9, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1F89, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1F02, 0x03B9, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1F8A, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1F03, 0x03B9, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1F8B, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1F04, 0x03B9, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1F8C, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1F05, 0x03B9, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1F8D, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1F06, 0x03B9, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1F8E, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1F07, 0x03B9, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1F8F, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1F20, 0x03B9, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1F90, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1F21, 0x03B9, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1F91, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1F22, 0x03B9, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1F92, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1F23, 0x03B9, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1F93, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1F24, 0x03B9, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1F94, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1F25, 0x03B9, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1F95, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1F26, 0x03B9, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1F96, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1F27, 0x03B9, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1F97, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1F20, 0x03B9, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1F98, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1F21, 0x03B9, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1F99, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1F22, 0x03B9, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1F9A, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1F23, 0x03B9, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1F9B, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1F24, 0x03B9, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1F9C, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1F25, 0x03B9, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1F9D, result.begin());
        BOOST_TEST(result == expected);
    }
}

// test_16)
{
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1F26, 0x03B9, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1F9E, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1F27, 0x03B9, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1F9F, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1F60, 0x03B9, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1FA0, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1F61, 0x03B9, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1FA1, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1F62, 0x03B9, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1FA2, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1F63, 0x03B9, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1FA3, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1F64, 0x03B9, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1FA4, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1F65, 0x03B9, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1FA5, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1F66, 0x03B9, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1FA6, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1F67, 0x03B9, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1FA7, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1F60, 0x03B9, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1FA8, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1F61, 0x03B9, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1FA9, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1F62, 0x03B9, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1FAA, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1F63, 0x03B9, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1FAB, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1F64, 0x03B9, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1FAC, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1F65, 0x03B9, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1FAD, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1F66, 0x03B9, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1FAE, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1F67, 0x03B9, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1FAF, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1F70, 0x03B9, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1FB2, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x03B1, 0x03B9, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1FB3, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x03AC, 0x03B9, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1FB4, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x03B1, 0x0342, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1FB6, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x03B1, 0x0342, 0x03B9, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1FB7, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1FB0, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1FB8, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1FB1, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1FB9, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1F70, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1FBA, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1F71, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1FBB, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x03B1, 0x03B9, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1FBC, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x03B9, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1FBE, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1F74, 0x03B9, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1FC2, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x03B7, 0x03B9, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1FC3, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x03AE, 0x03B9, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1FC4, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x03B7, 0x0342, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1FC6, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x03B7, 0x0342, 0x03B9, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1FC7, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1F72, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1FC8, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1F73, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1FC9, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1F74, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1FCA, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1F75, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1FCB, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x03B7, 0x03B9, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1FCC, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x03B9, 0x0308, 0x0300, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1FD2, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x03B9, 0x0308, 0x0301, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1FD3, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x03B9, 0x0342, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1FD6, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x03B9, 0x0308, 0x0342, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1FD7, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1FD0, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1FD8, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1FD1, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1FD9, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1F76, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1FDA, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1F77, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1FDB, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x03C5, 0x0308, 0x0300, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1FE2, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x03C5, 0x0308, 0x0301, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1FE3, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x03C1, 0x0313, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1FE4, result.begin());
        BOOST_TEST(result == expected);
    }
}

// test_17)
{
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x03C5, 0x0342, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1FE6, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x03C5, 0x0308, 0x0342, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1FE7, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1FE0, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1FE8, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1FE1, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1FE9, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1F7A, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1FEA, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1F7B, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1FEB, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1FE5, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1FEC, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1F7C, 0x03B9, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1FF2, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x03C9, 0x03B9, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1FF3, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x03CE, 0x03B9, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1FF4, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x03C9, 0x0342, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1FF6, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x03C9, 0x0342, 0x03B9, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1FF7, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1F78, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1FF8, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1F79, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1FF9, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1F7C, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1FFA, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1F7D, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1FFB, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x03C9, 0x03B9, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1FFC, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x03C9, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x2126, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x006B, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x212A, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x00E5, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x212B, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x214E, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x2132, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x2170, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x2160, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x2171, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x2161, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x2172, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x2162, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x2173, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x2163, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x2174, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x2164, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x2175, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x2165, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x2176, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x2166, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x2177, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x2167, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x2178, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x2168, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x2179, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x2169, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x217A, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x216A, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x217B, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x216B, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x217C, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x216C, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x217D, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x216D, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x217E, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x216E, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x217F, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x216F, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x2184, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x2183, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x24D0, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x24B6, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x24D1, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x24B7, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x24D2, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x24B8, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x24D3, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x24B9, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x24D4, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x24BA, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x24D5, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x24BB, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x24D6, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x24BC, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x24D7, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x24BD, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x24D8, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x24BE, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x24D9, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x24BF, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x24DA, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x24C0, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x24DB, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x24C1, result.begin());
        BOOST_TEST(result == expected);
    }
}

// test_18)
{
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x24DC, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x24C2, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x24DD, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x24C3, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x24DE, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x24C4, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x24DF, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x24C5, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x24E0, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x24C6, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x24E1, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x24C7, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x24E2, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x24C8, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x24E3, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x24C9, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x24E4, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x24CA, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x24E5, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x24CB, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x24E6, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x24CC, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x24E7, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x24CD, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x24E8, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x24CE, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x24E9, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x24CF, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x2C30, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x2C00, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x2C31, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x2C01, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x2C32, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x2C02, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x2C33, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x2C03, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x2C34, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x2C04, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x2C35, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x2C05, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x2C36, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x2C06, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x2C37, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x2C07, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x2C38, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x2C08, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x2C39, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x2C09, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x2C3A, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x2C0A, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x2C3B, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x2C0B, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x2C3C, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x2C0C, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x2C3D, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x2C0D, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x2C3E, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x2C0E, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x2C3F, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x2C0F, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x2C40, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x2C10, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x2C41, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x2C11, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x2C42, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x2C12, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x2C43, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x2C13, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x2C44, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x2C14, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x2C45, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x2C15, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x2C46, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x2C16, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x2C47, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x2C17, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x2C48, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x2C18, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x2C49, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x2C19, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x2C4A, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x2C1A, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x2C4B, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x2C1B, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x2C4C, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x2C1C, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x2C4D, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x2C1D, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x2C4E, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x2C1E, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x2C4F, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x2C1F, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x2C50, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x2C20, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x2C51, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x2C21, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x2C52, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x2C22, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x2C53, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x2C23, result.begin());
        BOOST_TEST(result == expected);
    }
}

// test_19)
{
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x2C54, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x2C24, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x2C55, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x2C25, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x2C56, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x2C26, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x2C57, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x2C27, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x2C58, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x2C28, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x2C59, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x2C29, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x2C5A, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x2C2A, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x2C5B, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x2C2B, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x2C5C, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x2C2C, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x2C5D, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x2C2D, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x2C5E, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x2C2E, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x2C5F, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x2C2F, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x2C61, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x2C60, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x026B, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x2C62, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1D7D, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x2C63, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x027D, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x2C64, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x2C68, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x2C67, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x2C6A, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x2C69, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x2C6C, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x2C6B, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x0251, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x2C6D, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x0271, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x2C6E, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x0250, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x2C6F, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x0252, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x2C70, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x2C73, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x2C72, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x2C76, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x2C75, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x023F, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x2C7E, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x0240, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x2C7F, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x2C81, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x2C80, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x2C83, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x2C82, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x2C85, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x2C84, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x2C87, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x2C86, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x2C89, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x2C88, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x2C8B, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x2C8A, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x2C8D, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x2C8C, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x2C8F, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x2C8E, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x2C91, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x2C90, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x2C93, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x2C92, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x2C95, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x2C94, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x2C97, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x2C96, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x2C99, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x2C98, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x2C9B, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x2C9A, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x2C9D, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x2C9C, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x2C9F, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x2C9E, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x2CA1, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x2CA0, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x2CA3, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x2CA2, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x2CA5, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x2CA4, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x2CA7, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x2CA6, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x2CA9, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x2CA8, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x2CAB, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x2CAA, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x2CAD, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x2CAC, result.begin());
        BOOST_TEST(result == expected);
    }
}

// test_20)
{
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x2CAF, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x2CAE, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x2CB1, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x2CB0, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x2CB3, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x2CB2, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x2CB5, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x2CB4, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x2CB7, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x2CB6, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x2CB9, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x2CB8, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x2CBB, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x2CBA, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x2CBD, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x2CBC, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x2CBF, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x2CBE, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x2CC1, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x2CC0, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x2CC3, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x2CC2, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x2CC5, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x2CC4, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x2CC7, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x2CC6, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x2CC9, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x2CC8, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x2CCB, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x2CCA, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x2CCD, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x2CCC, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x2CCF, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x2CCE, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x2CD1, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x2CD0, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x2CD3, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x2CD2, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x2CD5, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x2CD4, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x2CD7, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x2CD6, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x2CD9, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x2CD8, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x2CDB, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x2CDA, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x2CDD, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x2CDC, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x2CDF, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x2CDE, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x2CE1, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x2CE0, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x2CE3, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x2CE2, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x2CEC, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x2CEB, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x2CEE, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x2CED, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x2CF3, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x2CF2, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0xA641, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xA640, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0xA643, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xA642, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0xA645, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xA644, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0xA647, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xA646, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0xA649, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xA648, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0xA64B, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xA64A, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0xA64D, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xA64C, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0xA64F, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xA64E, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0xA651, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xA650, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0xA653, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xA652, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0xA655, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xA654, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0xA657, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xA656, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0xA659, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xA658, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0xA65B, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xA65A, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0xA65D, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xA65C, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0xA65F, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xA65E, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0xA661, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xA660, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0xA663, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xA662, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0xA665, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xA664, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0xA667, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xA666, result.begin());
        BOOST_TEST(result == expected);
    }
}

// test_21)
{
    {
        std::array<uint32_t, 3 + 1> const expected = { 0xA669, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xA668, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0xA66B, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xA66A, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0xA66D, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xA66C, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0xA681, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xA680, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0xA683, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xA682, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0xA685, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xA684, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0xA687, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xA686, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0xA689, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xA688, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0xA68B, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xA68A, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0xA68D, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xA68C, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0xA68F, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xA68E, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0xA691, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xA690, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0xA693, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xA692, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0xA695, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xA694, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0xA697, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xA696, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0xA699, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xA698, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0xA69B, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xA69A, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0xA723, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xA722, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0xA725, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xA724, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0xA727, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xA726, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0xA729, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xA728, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0xA72B, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xA72A, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0xA72D, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xA72C, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0xA72F, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xA72E, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0xA733, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xA732, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0xA735, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xA734, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0xA737, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xA736, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0xA739, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xA738, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0xA73B, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xA73A, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0xA73D, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xA73C, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0xA73F, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xA73E, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0xA741, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xA740, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0xA743, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xA742, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0xA745, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xA744, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0xA747, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xA746, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0xA749, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xA748, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0xA74B, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xA74A, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0xA74D, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xA74C, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0xA74F, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xA74E, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0xA751, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xA750, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0xA753, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xA752, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0xA755, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xA754, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0xA757, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xA756, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0xA759, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xA758, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0xA75B, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xA75A, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0xA75D, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xA75C, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0xA75F, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xA75E, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0xA761, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xA760, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0xA763, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xA762, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0xA765, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xA764, result.begin());
        BOOST_TEST(result == expected);
    }
}

// test_22)
{
    {
        std::array<uint32_t, 3 + 1> const expected = { 0xA767, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xA766, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0xA769, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xA768, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0xA76B, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xA76A, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0xA76D, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xA76C, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0xA76F, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xA76E, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0xA77A, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xA779, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0xA77C, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xA77B, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1D79, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xA77D, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0xA77F, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xA77E, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0xA781, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xA780, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0xA783, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xA782, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0xA785, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xA784, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0xA787, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xA786, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0xA78C, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xA78B, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x0265, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xA78D, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0xA791, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xA790, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0xA793, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xA792, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0xA797, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xA796, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0xA799, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xA798, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0xA79B, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xA79A, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0xA79D, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xA79C, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0xA79F, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xA79E, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0xA7A1, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xA7A0, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0xA7A3, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xA7A2, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0xA7A5, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xA7A4, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0xA7A7, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xA7A6, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0xA7A9, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xA7A8, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x0266, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xA7AA, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x025C, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xA7AB, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x0261, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xA7AC, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x026C, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xA7AD, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x026A, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xA7AE, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x029E, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xA7B0, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x0287, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xA7B1, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x029D, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xA7B2, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0xAB53, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xA7B3, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0xA7B5, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xA7B4, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0xA7B7, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xA7B6, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0xA7B9, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xA7B8, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0xA7BB, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xA7BA, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0xA7BD, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xA7BC, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0xA7BF, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xA7BE, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0xA7C1, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xA7C0, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0xA7C3, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xA7C2, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0xA794, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xA7C4, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x0282, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xA7C5, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1D8E, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xA7C6, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0xA7C8, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xA7C7, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0xA7CA, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xA7C9, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0xA7D1, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xA7D0, result.begin());
        BOOST_TEST(result == expected);
    }
}

// test_23)
{
    {
        std::array<uint32_t, 3 + 1> const expected = { 0xA7D7, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xA7D6, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0xA7D9, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xA7D8, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0xA7F6, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xA7F5, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x13A0, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xAB70, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x13A1, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xAB71, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x13A2, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xAB72, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x13A3, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xAB73, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x13A4, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xAB74, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x13A5, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xAB75, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x13A6, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xAB76, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x13A7, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xAB77, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x13A8, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xAB78, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x13A9, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xAB79, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x13AA, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xAB7A, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x13AB, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xAB7B, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x13AC, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xAB7C, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x13AD, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xAB7D, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x13AE, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xAB7E, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x13AF, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xAB7F, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x13B0, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xAB80, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x13B1, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xAB81, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x13B2, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xAB82, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x13B3, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xAB83, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x13B4, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xAB84, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x13B5, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xAB85, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x13B6, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xAB86, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x13B7, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xAB87, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x13B8, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xAB88, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x13B9, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xAB89, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x13BA, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xAB8A, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x13BB, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xAB8B, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x13BC, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xAB8C, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x13BD, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xAB8D, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x13BE, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xAB8E, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x13BF, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xAB8F, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x13C0, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xAB90, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x13C1, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xAB91, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x13C2, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xAB92, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x13C3, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xAB93, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x13C4, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xAB94, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x13C5, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xAB95, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x13C6, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xAB96, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x13C7, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xAB97, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x13C8, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xAB98, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x13C9, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xAB99, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x13CA, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xAB9A, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x13CB, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xAB9B, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x13CC, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xAB9C, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x13CD, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xAB9D, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x13CE, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xAB9E, result.begin());
        BOOST_TEST(result == expected);
    }
}

// test_24)
{
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x13CF, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xAB9F, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x13D0, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xABA0, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x13D1, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xABA1, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x13D2, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xABA2, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x13D3, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xABA3, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x13D4, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xABA4, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x13D5, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xABA5, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x13D6, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xABA6, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x13D7, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xABA7, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x13D8, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xABA8, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x13D9, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xABA9, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x13DA, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xABAA, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x13DB, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xABAB, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x13DC, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xABAC, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x13DD, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xABAD, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x13DE, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xABAE, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x13DF, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xABAF, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x13E0, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xABB0, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x13E1, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xABB1, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x13E2, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xABB2, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x13E3, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xABB3, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x13E4, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xABB4, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x13E5, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xABB5, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x13E6, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xABB6, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x13E7, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xABB7, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x13E8, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xABB8, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x13E9, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xABB9, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x13EA, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xABBA, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x13EB, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xABBB, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x13EC, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xABBC, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x13ED, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xABBD, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x13EE, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xABBE, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x13EF, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xABBF, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x0066, 0x0066, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xFB00, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x0066, 0x0069, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xFB01, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x0066, 0x006C, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xFB02, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x0066, 0x0066, 0x0069, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xFB03, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x0066, 0x0066, 0x006C, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xFB04, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x0073, 0x0074, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xFB05, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x0073, 0x0074, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xFB06, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x0574, 0x0576, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xFB13, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x0574, 0x0565, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xFB14, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x0574, 0x056B, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xFB15, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x057E, 0x0576, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xFB16, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x0574, 0x056D, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xFB17, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0xFF41, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xFF21, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0xFF42, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xFF22, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0xFF43, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xFF23, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0xFF44, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xFF24, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0xFF45, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xFF25, result.begin());
        BOOST_TEST(result == expected);
    }
}

// test_25)
{
    {
        std::array<uint32_t, 3 + 1> const expected = { 0xFF46, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xFF26, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0xFF47, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xFF27, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0xFF48, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xFF28, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0xFF49, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xFF29, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0xFF4A, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xFF2A, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0xFF4B, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xFF2B, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0xFF4C, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xFF2C, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0xFF4D, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xFF2D, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0xFF4E, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xFF2E, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0xFF4F, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xFF2F, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0xFF50, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xFF30, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0xFF51, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xFF31, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0xFF52, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xFF32, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0xFF53, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xFF33, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0xFF54, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xFF34, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0xFF55, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xFF35, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0xFF56, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xFF36, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0xFF57, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xFF37, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0xFF58, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xFF38, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0xFF59, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xFF39, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0xFF5A, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0xFF3A, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x10428, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x10400, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x10429, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x10401, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1042A, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x10402, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1042B, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x10403, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1042C, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x10404, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1042D, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x10405, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1042E, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x10406, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1042F, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x10407, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x10430, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x10408, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x10431, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x10409, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x10432, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1040A, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x10433, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1040B, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x10434, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1040C, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x10435, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1040D, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x10436, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1040E, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x10437, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1040F, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x10438, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x10410, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x10439, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x10411, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1043A, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x10412, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1043B, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x10413, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1043C, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x10414, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1043D, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x10415, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1043E, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x10416, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1043F, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x10417, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x10440, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x10418, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x10441, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x10419, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x10442, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1041A, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x10443, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1041B, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x10444, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1041C, result.begin());
        BOOST_TEST(result == expected);
    }
}

// test_26)
{
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x10445, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1041D, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x10446, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1041E, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x10447, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1041F, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x10448, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x10420, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x10449, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x10421, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1044A, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x10422, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1044B, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x10423, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1044C, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x10424, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1044D, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x10425, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1044E, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x10426, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1044F, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x10427, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x104D8, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x104B0, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x104D9, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x104B1, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x104DA, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x104B2, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x104DB, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x104B3, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x104DC, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x104B4, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x104DD, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x104B5, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x104DE, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x104B6, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x104DF, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x104B7, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x104E0, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x104B8, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x104E1, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x104B9, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x104E2, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x104BA, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x104E3, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x104BB, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x104E4, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x104BC, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x104E5, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x104BD, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x104E6, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x104BE, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x104E7, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x104BF, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x104E8, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x104C0, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x104E9, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x104C1, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x104EA, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x104C2, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x104EB, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x104C3, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x104EC, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x104C4, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x104ED, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x104C5, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x104EE, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x104C6, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x104EF, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x104C7, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x104F0, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x104C8, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x104F1, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x104C9, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x104F2, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x104CA, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x104F3, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x104CB, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x104F4, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x104CC, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x104F5, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x104CD, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x104F6, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x104CE, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x104F7, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x104CF, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x104F8, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x104D0, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x104F9, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x104D1, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x104FA, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x104D2, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x104FB, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x104D3, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x10597, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x10570, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x10598, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x10571, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x10599, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x10572, result.begin());
        BOOST_TEST(result == expected);
    }
}

// test_27)
{
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1059A, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x10573, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1059B, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x10574, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1059C, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x10575, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1059D, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x10576, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1059E, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x10577, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1059F, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x10578, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x105A0, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x10579, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x105A1, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1057A, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x105A3, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1057C, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x105A4, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1057D, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x105A5, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1057E, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x105A6, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1057F, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x105A7, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x10580, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x105A8, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x10581, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x105A9, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x10582, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x105AA, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x10583, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x105AB, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x10584, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x105AC, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x10585, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x105AD, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x10586, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x105AE, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x10587, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x105AF, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x10588, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x105B0, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x10589, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x105B1, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1058A, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x105B3, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1058C, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x105B4, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1058D, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x105B5, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1058E, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x105B6, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1058F, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x105B7, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x10590, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x105B8, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x10591, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x105B9, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x10592, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x105BB, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x10594, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x105BC, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x10595, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x10CC0, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x10C80, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x10CC1, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x10C81, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x10CC2, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x10C82, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x10CC3, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x10C83, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x10CC4, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x10C84, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x10CC5, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x10C85, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x10CC6, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x10C86, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x10CC7, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x10C87, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x10CC8, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x10C88, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x10CC9, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x10C89, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x10CCA, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x10C8A, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x10CCB, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x10C8B, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x10CCC, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x10C8C, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x10CCD, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x10C8D, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x10CCE, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x10C8E, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x10CCF, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x10C8F, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x10CD0, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x10C90, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x10CD1, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x10C91, result.begin());
        BOOST_TEST(result == expected);
    }
}

// test_28)
{
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x10CD2, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x10C92, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x10CD3, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x10C93, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x10CD4, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x10C94, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x10CD5, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x10C95, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x10CD6, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x10C96, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x10CD7, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x10C97, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x10CD8, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x10C98, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x10CD9, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x10C99, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x10CDA, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x10C9A, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x10CDB, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x10C9B, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x10CDC, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x10C9C, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x10CDD, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x10C9D, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x10CDE, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x10C9E, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x10CDF, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x10C9F, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x10CE0, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x10CA0, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x10CE1, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x10CA1, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x10CE2, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x10CA2, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x10CE3, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x10CA3, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x10CE4, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x10CA4, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x10CE5, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x10CA5, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x10CE6, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x10CA6, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x10CE7, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x10CA7, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x10CE8, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x10CA8, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x10CE9, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x10CA9, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x10CEA, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x10CAA, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x10CEB, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x10CAB, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x10CEC, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x10CAC, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x10CED, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x10CAD, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x10CEE, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x10CAE, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x10CEF, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x10CAF, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x10CF0, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x10CB0, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x10CF1, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x10CB1, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x10CF2, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x10CB2, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x118C0, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x118A0, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x118C1, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x118A1, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x118C2, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x118A2, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x118C3, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x118A3, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x118C4, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x118A4, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x118C5, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x118A5, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x118C6, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x118A6, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x118C7, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x118A7, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x118C8, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x118A8, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x118C9, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x118A9, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x118CA, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x118AA, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x118CB, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x118AB, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x118CC, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x118AC, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x118CD, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x118AD, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x118CE, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x118AE, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x118CF, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x118AF, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x118D0, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x118B0, result.begin());
        BOOST_TEST(result == expected);
    }
}

// test_29)
{
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x118D1, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x118B1, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x118D2, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x118B2, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x118D3, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x118B3, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x118D4, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x118B4, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x118D5, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x118B5, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x118D6, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x118B6, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x118D7, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x118B7, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x118D8, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x118B8, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x118D9, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x118B9, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x118DA, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x118BA, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x118DB, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x118BB, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x118DC, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x118BC, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x118DD, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x118BD, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x118DE, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x118BE, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x118DF, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x118BF, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x16E60, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x16E40, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x16E61, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x16E41, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x16E62, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x16E42, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x16E63, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x16E43, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x16E64, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x16E44, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x16E65, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x16E45, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x16E66, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x16E46, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x16E67, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x16E47, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x16E68, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x16E48, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x16E69, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x16E49, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x16E6A, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x16E4A, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x16E6B, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x16E4B, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x16E6C, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x16E4C, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x16E6D, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x16E4D, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x16E6E, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x16E4E, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x16E6F, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x16E4F, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x16E70, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x16E50, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x16E71, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x16E51, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x16E72, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x16E52, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x16E73, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x16E53, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x16E74, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x16E54, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x16E75, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x16E55, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x16E76, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x16E56, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x16E77, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x16E57, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x16E78, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x16E58, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x16E79, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x16E59, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x16E7A, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x16E5A, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x16E7B, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x16E5B, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x16E7C, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x16E5C, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x16E7D, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x16E5D, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x16E7E, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x16E5E, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x16E7F, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x16E5F, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1E922, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1E900, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1E923, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1E901, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1E924, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1E902, result.begin());
        BOOST_TEST(result == expected);
    }
}

// test_30)
{
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1E925, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1E903, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1E926, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1E904, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1E927, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1E905, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1E928, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1E906, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1E929, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1E907, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1E92A, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1E908, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1E92B, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1E909, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1E92C, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1E90A, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1E92D, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1E90B, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1E92E, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1E90C, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1E92F, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1E90D, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1E930, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1E90E, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1E931, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1E90F, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1E932, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1E910, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1E933, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1E911, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1E934, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1E912, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1E935, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1E913, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1E936, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1E914, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1E937, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1E915, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1E938, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1E916, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1E939, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1E917, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1E93A, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1E918, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1E93B, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1E919, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1E93C, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1E91A, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1E93D, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1E91B, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1E93E, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1E91C, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1E93F, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1E91D, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1E940, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1E91E, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1E941, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1E91F, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1E942, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1E920, result.begin());
        BOOST_TEST(result == expected);
    }
    {
        std::array<uint32_t, 3 + 1> const expected = { 0x1E943, 0, 0, 0 };
        std::array<uint32_t, 3 + 1> result = { 0, 0, 0, 0 };
        boost::parser::detail::case_fold(0x1E921, result.begin());
        BOOST_TEST(result == expected);
    }
}

// misses
{
    char32_t next_cp = 0;
    for (char32_t const * it = cps; it != std::end(cps); ++it) {
        for (char32_t cp = next_cp; cp < *it; ++cp) {
            std::array<uint32_t, 3 + 1> const expected = { cp, 0 };
            std::array<uint32_t, 3 + 1> result = { 0 };
            auto const first = result.data();
            auto const last = boost::parser::detail::case_fold(cp, first);
            BOOST_TEST(std::equal(first, last, expected.begin()));
            BOOST_TEST(result == expected);
        }
        next_cp = *it + 1;
    }
}

return boost::report_errors();
}
