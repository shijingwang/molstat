#define max_atoms 1024
#define max_bonds 1024
#define fpf_boolean 3001
#define fpf_decimal 3002
#define linesize 500 //give a large enough number to read molfile. 
#define  max_matchpath_length 256
#define  pmCheckMol   1001
#define  pmMatchMol   1002
#define  tweaklevel   2

#define  rs_sar       2001   // ring search mode: SAR = set of all rings
#define  rs_ssr       2002   //                   SSR = set of small rings

#define molfileend      "$$$$"
#define molend          "M  END"
#define  max_ringsize   128
#define  max_rings      1024
#define  max_fg         256 // must be a multiple of 32
#define  used_fg        205 // v0.5; number of functional groups currently used (204 in v0.5)
#define  max_neighbors  20  // new in v0.2h (16); increased in v0.4d (20)
#define  TAB            #26
#define  max_matchpath_length  256

#define  btopo_any         0  // bond topology, new in v0.3d
#define  btopo_ring        1  //
#define  btopo_chain       2  //
#define  btopo_always_any  3  // even in "strict mode"
#define  btopo_excess_rc   4  // bond in candidate must be included in _more_ rings than
                         // the matching bond in the query ==> specific search for
                         // annulated systems
#define  btopo_exact_rc    5  // bond in query and candidate must have same ring count
#define  bstereo_any       0 // new in v0.3d, any E/Z isomer matches (for double bonds)
#define  bstereo_xyz       1  // E/Z match is checked by using XYZ coordinates of the atoms
#define  bstereo_up        11 // new in v0.3f, flags for single bonds
#define  bstereo_down      16 //

#define  fp_blocksize      62   // v0.3m; 1 bit may be used for sign (e.g. by PHP), 1 bit for "exact hit"
#define  ssr_vringsize     12  // v0.3n; max. ring size in SSR mode
#define  max_recursion_depth        500000  // v0.3p
#define  max_match_recursion_depth  500000  // v0.3p
#define  max_ndl_gmmsize   128  // v0.4a
#define  max_hay_gmmsize   1024 // v0.4a
#define  pi                3.1415926535897932384626433832795028841971693993751

#define  max_atomicnum            118 // v0.4
#define  max_fragpath_length      8   // v0.4
#define  min_fragpath_length      3   // v0.4
#define  max_ringfragpath_length  12  // should be identical with ssr_vringsize
#define  max_fragstr_length       64  // v0.4; must be at least 4 * max_ringfragpath_length
#define  hfpsize                  512 // v0.4
#define  n_hfpbits                2   // v0.4; currently supported: 2 or 3
#define  max_n_matches            64  // v0.4b; maximum number of matches of one and the same needle+haystack pair
#define  lang_en                  1   // v0.5
//  lang_de                 = 2;   // v0.5
#define  max_fgpos                127  // v0.5 max. number of located funtional groups per fg type

  // Definitions for functional groups:
#define   fg_cation                            1
#define   fg_anion                             2
#define   fg_carbonyl                          3
#define   fg_aldehyde                          4
#define   fg_ketone                            5
#define   fg_thiocarbonyl                      6
#define   fg_thioaldehyde                      7
#define   fg_thioketone                        8
#define   fg_imine                             9
#define   fg_hydrazone                        10
#define   fg_semicarbazone                    11
#define   fg_thiosemicarbazone                12
#define   fg_oxime                            13
#define   fg_oxime_ether                      14
#define   fg_ketene                           15
#define   fg_ketene_acetal_deriv              16
#define   fg_carbonyl_hydrate                 17
#define   fg_hemiacetal                       18
#define   fg_acetal                           19
#define   fg_hemiaminal                       20
#define   fg_aminal                           21
#define   fg_thiohemiaminal                   22
#define   fg_thioacetal                       23
#define   fg_enamine                          24
#define   fg_enol                             25
#define   fg_enolether                        26
#define   fg_hydroxy                          27
#define   fg_alcohol                          28
#define   fg_prim_alcohol                     29
#define   fg_sec_alcohol                      30
#define   fg_tert_alcohol                     31
#define   fg_1_2_diol                         32
#define   fg_1_2_aminoalcohol                 33
#define   fg_phenol                           34
#define   fg_1_2_diphenol                    35
#define   fg_enediol                         36
#define   fg_ether                           37
#define   fg_dialkylether                    38
#define   fg_alkylarylether                  39
#define   fg_diarylether                     40
#define   fg_thioether                       41
#define   fg_disulfide                       42
#define   fg_peroxide                        43
#define   fg_hydroperoxide                   44
#define   fg_hydrazine                       45
#define   fg_hydroxylamine                   46
#define   fg_amine                           47
#define   fg_prim_amine                      48
#define   fg_prim_aliph_amine                49
#define   fg_prim_arom_amine                 50
#define   fg_sec_amine                       51
#define   fg_sec_aliph_amine                 52
#define   fg_sec_mixed_amine                 53
#define   fg_sec_arom_amine                  54
#define   fg_tert_amine                      55
#define   fg_tert_aliph_amine                56
#define   fg_tert_mixed_amine                57
#define   fg_tert_arom_amine                 58
#define   fg_quart_ammonium                  59
#define   fg_n_oxide                         60
#define   fg_halogen_deriv                   61
#define   fg_alkyl_halide                    62
#define   fg_alkyl_fluoride                  63
#define   fg_alkyl_chloride                  64
#define   fg_alkyl_bromide                   65
#define   fg_alkyl_iodide                    66
#define   fg_aryl_halide                     67
#define   fg_aryl_fluoride                   68
#define   fg_aryl_chloride                   69
#define   fg_aryl_bromide                    70
#define   fg_aryl_iodide                     71
#define   fg_organometallic                  72
#define   fg_organolithium                   73
#define   fg_organomagnesium                 74
#define   fg_carboxylic_acid_deriv           75
#define   fg_carboxylic_acid                 76
#define   fg_carboxylic_acid_salt            77
#define   fg_carboxylic_acid_ester           78
#define   fg_lactone                         79
#define   fg_carboxylic_acid_amide           80
#define   fg_carboxylic_acid_prim_amide      81
#define   fg_carboxylic_acid_sec_amide       82
#define   fg_carboxylic_acid_tert_amide      83
#define   fg_lactam                          84
#define   fg_carboxylic_acid_hydrazide       85
#define   fg_carboxylic_acid_azide           86
#define   fg_hydroxamic_acid                 87
#define   fg_carboxylic_acid_amidine         88
#define   fg_carboxylic_acid_amidrazone      89
#define   fg_nitrile                         90
#define   fg_acyl_halide                     91
#define   fg_acyl_fluoride                   92
#define   fg_acyl_chloride                   93
#define   fg_acyl_bromide                    94
#define   fg_acyl_iodide                     95
#define   fg_acyl_cyanide                    96
#define   fg_imido_ester                     97
#define   fg_imidoyl_halide                  98
#define   fg_thiocarboxylic_acid_deriv       99
#define   fg_thiocarboxylic_acid            100
#define   fg_thiocarboxylic_acid_ester      101
#define   fg_thiolactone                    102
#define   fg_thiocarboxylic_acid_amide      103
#define   fg_thiolactam                     104
#define   fg_imido_thioester                105
#define   fg_oxohetarene                    106
#define   fg_thioxohetarene                 107
#define   fg_iminohetarene                  108
#define   fg_orthocarboxylic_acid_deriv     109
#define   fg_carboxylic_acid_orthoester     110
#define   fg_carboxylic_acid_amide_acetal   111
#define   fg_carboxylic_acid_anhydride      112
#define   fg_carboxylic_acid_imide          113
#define   fg_carboxylic_acid_unsubst_imide  114
#define   fg_carboxylic_acid_subst_imide    115
#define   fg_co2_deriv                      116
#define   fg_carbonic_acid_deriv            117
#define   fg_carbonic_acid_monoester        118
#define   fg_carbonic_acid_diester          119
#define   fg_carbonic_acid_ester_halide     120
#define   fg_thiocarbonic_acid_deriv        121
#define   fg_thiocarbonic_acid_monoester    122
#define   fg_thiocarbonic_acid_diester      123
#define   fg_thiocarbonic_acid_ester_halide 124
#define   fg_carbamic_acid_deriv            125
#define   fg_carbamic_acid                  126
#define   fg_carbamic_acid_ester            127
#define   fg_carbamic_acid_halide           128
#define   fg_thiocarbamic_acid_deriv        129
#define   fg_thiocarbamic_acid              130
#define   fg_thiocarbamic_acid_ester        131
#define   fg_thiocarbamic_acid_halide       132
#define   fg_urea                            133
#define   fg_isourea                         134
#define   fg_thiourea                        135
#define   fg_isothiourea                     136
#define   fg_guanidine                       137
#define   fg_semicarbazide                   138
#define   fg_thiosemicarbazide               139
#define   fg_azide                           140
#define   fg_azo_compound                    141
#define   fg_diazonium_salt                  142
#define   fg_isonitrile                      143
#define   fg_cyanate                         144
#define   fg_isocyanate                      145
#define   fg_thiocyanate                     146
#define   fg_isothiocyanate                  147
#define   fg_carbodiimide                    148
#define   fg_nitroso_compound                149
#define   fg_nitro_compound                  150
#define   fg_nitrite                         151
#define   fg_nitrate                         152
#define   fg_sulfuric_acid_deriv             153
#define   fg_sulfuric_acid                   154
#define   fg_sulfuric_acid_monoester         155
#define   fg_sulfuric_acid_diester           156
#define   fg_sulfuric_acid_amide_ester       157
#define   fg_sulfuric_acid_amide             158
#define   fg_sulfuric_acid_diamide           159
#define   fg_sulfuryl_halide                 160
#define   fg_sulfonic_acid_deriv             161
#define   fg_sulfonic_acid                   162
#define   fg_sulfonic_acid_ester             163
#define   fg_sulfonamide                     164
#define   fg_sulfonyl_halide                 165
#define   fg_sulfone                         166
#define   fg_sulfoxide                       167
#define   fg_sulfinic_acid_deriv             168
#define   fg_sulfinic_acid                   169
#define   fg_sulfinic_acid_ester             170
#define   fg_sulfinic_acid_halide            171
#define   fg_sulfinic_acid_amide             172
#define   fg_sulfenic_acid_deriv             173
#define   fg_sulfenic_acid                   174
#define   fg_sulfenic_acid_ester             175
#define   fg_sulfenic_acid_halide            176
#define   fg_sulfenic_acid_amide             177
#define   fg_thiol                           178
#define   fg_alkylthiol                      179
#define   fg_arylthiol                       180
#define   fg_phosphoric_acid_deriv           181
#define   fg_phosphoric_acid                 182
#define   fg_phosphoric_acid_ester           183
#define   fg_phosphoric_acid_halide          184
#define   fg_phosphoric_acid_amide           185
#define   fg_thiophosphoric_acid_deriv       186
#define   fg_thiophosphoric_acid             187
#define   fg_thiophosphoric_acid_ester       188
#define   fg_thiophosphoric_acid_halide      189
#define   fg_thiophosphoric_acid_amide       190
#define   fg_phosphonic_acid_deriv           191
#define   fg_phosphonic_acid                 192
#define   fg_phosphonic_acid_ester           193
#define   fg_phosphine                       194
#define   fg_phosphinoxide                   195
#define   fg_boronic_acid_deriv              196
#define   fg_boronic_acid                    197
#define   fg_boronic_acid_ester              198
#define   fg_alkene                          199
#define   fg_alkyne                          200
#define   fg_aromatic                        201
#define   fg_heterocycle                     202
#define   fg_alpha_aminoacid                 203
#define   fg_alpha_hydroxyacid               204

