/* Copyright 2019 Axel Huebl, Maxence Thevenet
 *
 * This file is part of WarpX.
 *
 * License: BSD-3-Clause-LBNL
 */
// This script was automatically generated!
// Edit dev/Source/Utils/write_atomic_data_cpp.py instead!
#ifndef WARPX_IONIZATION_TABLE_H_
#define WARPX_IONIZATION_TABLE_H_

#include <AMReX_AmrCore.H>
#include <AMReX_REAL.H>
#include <map>

std::map<std::string, int> ion_map_ids = {
    {"H", 0},
    {"He", 1},
    {"Li", 2},
    {"Be", 3},
    {"B", 4},
    {"C", 5},
    {"N", 6},
    {"O", 7},
    {"F", 8},
    {"Ne", 9},
    {"Na", 10},
    {"Mg", 11},
    {"Al", 12},
    {"Si", 13},
    {"P", 14},
    {"S", 15},
    {"Cl", 16},
    {"Ar", 17},
    {"Kr", 18},
    {"Rb", 19},
    {"Xe", 20},
    {"Rn", 21} };

const int nelements = 22;

const int ion_atomic_numbers[nelements] = {
    1, 2, 3, 4, 5, 6, 7, 8, 9, 10,
    11, 12, 13, 14, 15, 16, 17, 18, 36, 37,
    54, 86};

const int ion_energy_offsets[nelements] = {
    0, 1, 3, 6, 10, 15, 21, 28, 36, 45,
    55, 66, 78, 91, 105, 120, 136, 153, 171, 207,
    244, 298};

const int energies_tab_length = 384;

const amrex::Real table_ionization_energies[energies_tab_length]{
    // H
    amrex::Real(13.59843449),
    // He
    amrex::Real(24.58738880), amrex::Real(54.4177650),
    // Li
    amrex::Real(5.39171495), amrex::Real(75.6400964), amrex::Real(122.4543581),
    // Be
    amrex::Real(9.322699), amrex::Real(18.21115), amrex::Real(153.896203),
    amrex::Real(217.7185843),
    // B
    amrex::Real(8.298019), amrex::Real(25.15483), amrex::Real(37.93058),
    amrex::Real(259.3715), amrex::Real(340.226020),
    // C
    amrex::Real(11.2602880), amrex::Real(24.383154), amrex::Real(47.88778),
    amrex::Real(64.49352), amrex::Real(392.090515), amrex::Real(489.993194),
    // N
    amrex::Real(14.53413), amrex::Real(29.60125), amrex::Real(47.4453),
    amrex::Real(77.4735), amrex::Real(97.8901), amrex::Real(552.06732),
    amrex::Real(667.046116),
    // O
    amrex::Real(13.618055), amrex::Real(35.12112), amrex::Real(54.93554),
    amrex::Real(77.41350), amrex::Real(113.8990), amrex::Real(138.1189),
    amrex::Real(739.32682), amrex::Real(871.40988),
    // F
    amrex::Real(17.42282), amrex::Real(34.97081), amrex::Real(62.70798),
    amrex::Real(87.175), amrex::Real(114.249), amrex::Real(157.16311),
    amrex::Real(185.1868), amrex::Real(953.89804), amrex::Real(1103.11747),
    // Ne
    amrex::Real(21.564540), amrex::Real(40.96297), amrex::Real(63.4233),
    amrex::Real(97.1900), amrex::Real(126.247), amrex::Real(157.934),
    amrex::Real(207.271), amrex::Real(239.0970), amrex::Real(1195.80783),
    amrex::Real(1362.19915),
    // Na
    amrex::Real(5.1390769), amrex::Real(47.28636), amrex::Real(71.6200),
    amrex::Real(98.936), amrex::Real(138.404), amrex::Real(172.23),
    amrex::Real(208.504), amrex::Real(264.192), amrex::Real(299.856),
    amrex::Real(1465.13449), amrex::Real(1648.70218),
    // Mg
    amrex::Real(7.646236), amrex::Real(15.035271), amrex::Real(80.1436),
    amrex::Real(109.2654), amrex::Real(141.33), amrex::Real(186.76),
    amrex::Real(225.02), amrex::Real(265.924), amrex::Real(327.99),
    amrex::Real(367.489), amrex::Real(1761.80487), amrex::Real(1962.66365),
    // Al
    amrex::Real(5.985769), amrex::Real(18.82855), amrex::Real(28.447642),
    amrex::Real(119.9924), amrex::Real(153.8252), amrex::Real(190.49),
    amrex::Real(241.76), amrex::Real(284.64), amrex::Real(330.21),
    amrex::Real(398.65), amrex::Real(442.005), amrex::Real(2085.97700),
    amrex::Real(2304.14005),
    // Si
    amrex::Real(8.15168), amrex::Real(16.34585), amrex::Real(33.49300),
    amrex::Real(45.14179), amrex::Real(166.767), amrex::Real(205.279),
    amrex::Real(246.57), amrex::Real(303.59), amrex::Real(351.28),
    amrex::Real(401.38), amrex::Real(476.273), amrex::Real(523.415),
    amrex::Real(2437.65813), amrex::Real(2673.17753),
    // P
    amrex::Real(10.486686), amrex::Real(19.76949), amrex::Real(30.20264),
    amrex::Real(51.44387), amrex::Real(65.02511), amrex::Real(220.430),
    amrex::Real(263.57), amrex::Real(309.60), amrex::Real(372.31),
    amrex::Real(424.40), amrex::Real(479.44), amrex::Real(560.62),
    amrex::Real(611.741), amrex::Real(2816.90876), amrex::Real(3069.8415),
    // S
    amrex::Real(10.36001), amrex::Real(23.33788), amrex::Real(34.86),
    amrex::Real(47.222), amrex::Real(72.5945), amrex::Real(88.0529),
    amrex::Real(280.954), amrex::Real(328.794), amrex::Real(379.84),
    amrex::Real(447.7), amrex::Real(504.55), amrex::Real(564.41),
    amrex::Real(651.96), amrex::Real(706.994), amrex::Real(3223.7807),
    amrex::Real(3494.1879),
    // Cl
    amrex::Real(12.967632), amrex::Real(23.81364), amrex::Real(39.80),
    amrex::Real(53.24), amrex::Real(67.68), amrex::Real(96.94),
    amrex::Real(114.2013), amrex::Real(348.306), amrex::Real(400.851),
    amrex::Real(456.7), amrex::Real(530.0), amrex::Real(591.58),
    amrex::Real(656.30), amrex::Real(750.23), amrex::Real(809.198),
    amrex::Real(3658.3437), amrex::Real(3946.2909),
    // Ar
    amrex::Real(15.7596117), amrex::Real(27.62967), amrex::Real(40.735),
    amrex::Real(59.58), amrex::Real(74.84), amrex::Real(91.290),
    amrex::Real(124.41), amrex::Real(143.4567), amrex::Real(422.60),
    amrex::Real(479.76), amrex::Real(540.4), amrex::Real(619.0),
    amrex::Real(685.5), amrex::Real(755.13), amrex::Real(855.5),
    amrex::Real(918.375), amrex::Real(4120.6656), amrex::Real(4426.2228),
    // Kr
    amrex::Real(13.9996053), amrex::Real(24.35984), amrex::Real(35.838),
    amrex::Real(50.85), amrex::Real(64.69), amrex::Real(78.49),
    amrex::Real(109.13), amrex::Real(125.802), amrex::Real(233.0),
    amrex::Real(268), amrex::Real(308), amrex::Real(350), amrex::Real(391),
    amrex::Real(446), amrex::Real(492), amrex::Real(540), amrex::Real(591),
    amrex::Real(640), amrex::Real(785), amrex::Real(831.6), amrex::Real(882.8),
    amrex::Real(945), amrex::Real(999.0), amrex::Real(1042),
    amrex::Real(1155.0), amrex::Real(1205.23), amrex::Real(2928.9),
    amrex::Real(3072), amrex::Real(3228), amrex::Real(3380), amrex::Real(3584),
    amrex::Real(3752.0), amrex::Real(3971), amrex::Real(4109.083),
    amrex::Real(17296.420), amrex::Real(17936.209),
    // Rb
    amrex::Real(4.1771280), amrex::Real(27.28954), amrex::Real(39.247),
    amrex::Real(52.20), amrex::Real(68.44), amrex::Real(82.9),
    amrex::Real(98.67), amrex::Real(132.79), amrex::Real(150.628),
    amrex::Real(277.12), amrex::Real(313.1), amrex::Real(356.0),
    amrex::Real(400), amrex::Real(443), amrex::Real(502), amrex::Real(550),
    amrex::Real(601), amrex::Real(654), amrex::Real(706.0), amrex::Real(857),
    amrex::Real(905.3), amrex::Real(958.9), amrex::Real(1024),
    amrex::Real(1080), amrex::Real(1125), amrex::Real(1242.5),
    amrex::Real(1294.57), amrex::Real(3133.3), amrex::Real(3281),
    amrex::Real(3443), amrex::Real(3600), amrex::Real(3815), amrex::Real(3988),
    amrex::Real(4214), amrex::Real(4356.865), amrex::Real(18305.884),
    amrex::Real(18965.516),
    // Xe
    amrex::Real(12.1298436), amrex::Real(20.975), amrex::Real(31.05),
    amrex::Real(42.20), amrex::Real(54.1), amrex::Real(66.703),
    amrex::Real(91.6), amrex::Real(105.9778), amrex::Real(179.84),
    amrex::Real(202.0), amrex::Real(229.02), amrex::Real(255.0),
    amrex::Real(281), amrex::Real(314), amrex::Real(343), amrex::Real(374),
    amrex::Real(404), amrex::Real(434), amrex::Real(549), amrex::Real(582),
    amrex::Real(616), amrex::Real(650), amrex::Real(700), amrex::Real(736),
    amrex::Real(818), amrex::Real(857.0), amrex::Real(1493), amrex::Real(1571),
    amrex::Real(1653), amrex::Real(1742), amrex::Real(1826), amrex::Real(1919),
    amrex::Real(2023), amrex::Real(2113), amrex::Real(2209), amrex::Real(2300),
    amrex::Real(2556), amrex::Real(2637), amrex::Real(2726), amrex::Real(2811),
    amrex::Real(2975), amrex::Real(3068), amrex::Real(3243),
    amrex::Real(3333.8), amrex::Real(7660), amrex::Real(7889),
    amrex::Real(8144), amrex::Real(8382), amrex::Real(8971),
    amrex::Real(9243), amrex::Real(9581), amrex::Real(9810.37),
    amrex::Real(40271.724), amrex::Real(41299.71),
    // Rn
    amrex::Real(10.74850), amrex::Real(21.4), amrex::Real(29.4),
    amrex::Real(36.9), amrex::Real(52.9), amrex::Real(64.0), amrex::Real(88.0),
    amrex::Real(102.0), amrex::Real(154.0), amrex::Real(173.9),
    amrex::Real(195.0), amrex::Real(218.0), amrex::Real(240), amrex::Real(264),
    amrex::Real(293), amrex::Real(317), amrex::Real(342), amrex::Real(367),
    amrex::Real(488), amrex::Real(520), amrex::Real(550), amrex::Real(580),
    amrex::Real(640), amrex::Real(680), amrex::Real(760), amrex::Real(800),
    amrex::Real(850), amrex::Real(920), amrex::Real(980), amrex::Real(1050),
    amrex::Real(1110), amrex::Real(1180), amrex::Real(1250), amrex::Real(1310),
    amrex::Real(1390), amrex::Real(1460), amrex::Real(1520), amrex::Real(1590),
    amrex::Real(1660), amrex::Real(1720), amrex::Real(2033), amrex::Real(2094),
    amrex::Real(2158), amrex::Real(2227), amrex::Real(2293), amrex::Real(2357),
    amrex::Real(2467), amrex::Real(2535), amrex::Real(2606), amrex::Real(2674),
    amrex::Real(2944), amrex::Real(3010), amrex::Real(3082), amrex::Real(3149),
    amrex::Real(3433), amrex::Real(3510), amrex::Real(3699), amrex::Real(3777),
    amrex::Real(6169), amrex::Real(6318), amrex::Real(6476), amrex::Real(6646),
    amrex::Real(6807), amrex::Real(6964), amrex::Real(7283), amrex::Real(7450),
    amrex::Real(7630), amrex::Real(7800), amrex::Real(8260), amrex::Real(8410),
    amrex::Real(8570), amrex::Real(8710), amrex::Real(9610), amrex::Real(9780),
    amrex::Real(10120), amrex::Real(10290), amrex::Real(21770),
    amrex::Real(22160), amrex::Real(22600), amrex::Real(22990),
    amrex::Real(26310), amrex::Real(26830), amrex::Real(27490),
    amrex::Real(27903.1), amrex::Real(110842.0), amrex::Real(112843.7)
};

#endif // #ifndef WARPX_IONIZATION_TABLE_H_
