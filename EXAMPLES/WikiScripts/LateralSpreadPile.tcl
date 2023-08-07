################################################################
#                                                              #
# p-y spring and beam model for lateral spreading analysis the #
# p-y curve parameters have been reduced to account for the    #
# presence of the weaker liquefied layer                       #
#                                                              #
# Pile parameters:  D = 0.6, E =3e+007, nu =   0.3             #
# Soil parameters:  z =  10, T =   2, phi =  36, gamma =  17   #
#                                                              #
# Created By:  Chris McGann                                    #
#              Pedro Arduino                                   #
#            --University of Washington--                      #
#                                                              #
# ---> Basic units are kN and m                                #
#                                                              #
################################################################

wipe
set outDir ./Output
file mkdir $outDir

#-----------------------------------------------------------------------------------------
#  1. CREATE P-Y SPRING NODES AND BOUNDARY CONDITIONS
#-----------------------------------------------------------------------------------------
model BasicBuilder -ndm 3 -ndf 3

# nodes for soil end of the spring
node  2  0.000 0.000 39.50000
node  3  0.000 0.000 39.00000
node  4  0.000 0.000 38.50000
node  5  0.000 0.000 38.00000
node  6  0.000 0.000 37.50000
node  7  0.000 0.000 37.00000
node  8  0.000 0.000 36.50000
node  9  0.000 0.000 36.00000
node 10  0.000 0.000 35.50000
node 11  0.000 0.000 35.00000
node 12  0.000 0.000 34.50000
node 13  0.000 0.000 34.00000
node 14  0.000 0.000 33.50000
node 15  0.000 0.000 33.00000
node 16  0.000 0.000 32.50000
node 17  0.000 0.000 32.00000
node 18  0.000 0.000 31.50000
node 19  0.000 0.000 31.00000
node 20  0.000 0.000 30.50000
node 21  0.000 0.000 30.00000
node 22  0.000 0.000 29.50000
node 23  0.000 0.000 29.00000
node 24  0.000 0.000 28.50000
node 25  0.000 0.000 28.00000
node 26  0.000 0.000 27.50000
node 27  0.000 0.000 27.00000
node 28  0.000 0.000 26.50000
node 29  0.000 0.000 26.00000
node 30  0.000 0.000 25.50000
node 31  0.000 0.000 25.00000
node 32  0.000 0.000 24.50000
node 33  0.000 0.000 24.00000
node 34  0.000 0.000 23.50000
node 35  0.000 0.000 23.00000
node 36  0.000 0.000 22.50000
node 37  0.000 0.000 22.00000
node 38  0.000 0.000 21.50000
node 39  0.000 0.000 21.00000
node 40  0.000 0.000 20.50000
node 41  0.000 0.000 20.00000
node 42  0.000 0.000 19.50000
node 43  0.000 0.000 19.00000
node 44  0.000 0.000 18.50000
node 45  0.000 0.000 18.00000
node 46  0.000 0.000 17.50000
node 47  0.000 0.000 17.00000
node 48  0.000 0.000 16.50000
node 49  0.000 0.000 16.00000
node 50  0.000 0.000 15.50000
node 51  0.000 0.000 15.00000
node 52  0.000 0.000 14.50000
node 53  0.000 0.000 14.00000
node 54  0.000 0.000 13.50000
node 55  0.000 0.000 13.00000
node 56  0.000 0.000 12.50000
node 57  0.000 0.000 12.00000
node 58  0.000 0.000 11.50000
node 59  0.000 0.000 11.00000
node 60  0.000 0.000 10.50000
node 61  0.000 0.000 10.00000
node 62  0.000 0.000 9.50000
node 63  0.000 0.000 9.00000
node 64  0.000 0.000 8.50000
node 65  0.000 0.000 8.00000
node 66  0.000 0.000 7.50000
node 67  0.000 0.000 7.00000
node 68  0.000 0.000 6.50000
node 69  0.000 0.000 6.00000
node 70  0.000 0.000 5.50000
node 71  0.000 0.000 5.00000
node 72  0.000 0.000 4.50000
node 73  0.000 0.000 4.00000
node 74  0.000 0.000 3.50000
node 75  0.000 0.000 3.00000
node 76  0.000 0.000 2.50000
node 77  0.000 0.000 2.00000
node 78  0.000 0.000 1.50000
node 79  0.000 0.000 1.00000
node 80  0.000 0.000 0.50000
node 81  0.000 0.000 0.00000

# nodes for pile end of the spring
node 202  0.000 0.000 39.50000
node 203  0.000 0.000 39.00000
node 204  0.000 0.000 38.50000
node 205  0.000 0.000 38.00000
node 206  0.000 0.000 37.50000
node 207  0.000 0.000 37.00000
node 208  0.000 0.000 36.50000
node 209  0.000 0.000 36.00000
node 210  0.000 0.000 35.50000
node 211  0.000 0.000 35.00000
node 212  0.000 0.000 34.50000
node 213  0.000 0.000 34.00000
node 214  0.000 0.000 33.50000
node 215  0.000 0.000 33.00000
node 216  0.000 0.000 32.50000
node 217  0.000 0.000 32.00000
node 218  0.000 0.000 31.50000
node 219  0.000 0.000 31.00000
node 220  0.000 0.000 30.50000
node 221  0.000 0.000 30.00000
node 222  0.000 0.000 29.50000
node 223  0.000 0.000 29.00000
node 224  0.000 0.000 28.50000
node 225  0.000 0.000 28.00000
node 226  0.000 0.000 27.50000
node 227  0.000 0.000 27.00000
node 228  0.000 0.000 26.50000
node 229  0.000 0.000 26.00000
node 230  0.000 0.000 25.50000
node 231  0.000 0.000 25.00000
node 232  0.000 0.000 24.50000
node 233  0.000 0.000 24.00000
node 234  0.000 0.000 23.50000
node 235  0.000 0.000 23.00000
node 236  0.000 0.000 22.50000
node 237  0.000 0.000 22.00000
node 238  0.000 0.000 21.50000
node 239  0.000 0.000 21.00000
node 240  0.000 0.000 20.50000
node 241  0.000 0.000 20.00000
node 242  0.000 0.000 19.50000
node 243  0.000 0.000 19.00000
node 244  0.000 0.000 18.50000
node 245  0.000 0.000 18.00000
node 246  0.000 0.000 17.50000
node 247  0.000 0.000 17.00000
node 248  0.000 0.000 16.50000
node 249  0.000 0.000 16.00000
node 250  0.000 0.000 15.50000
node 251  0.000 0.000 15.00000
node 252  0.000 0.000 14.50000
node 253  0.000 0.000 14.00000
node 254  0.000 0.000 13.50000
node 255  0.000 0.000 13.00000
node 256  0.000 0.000 12.50000
node 257  0.000 0.000 12.00000
node 258  0.000 0.000 11.50000
node 259  0.000 0.000 11.00000
node 260  0.000 0.000 10.50000
node 261  0.000 0.000 10.00000
node 262  0.000 0.000 9.50000
node 263  0.000 0.000 9.00000
node 264  0.000 0.000 8.50000
node 265  0.000 0.000 8.00000
node 266  0.000 0.000 7.50000
node 267  0.000 0.000 7.00000
node 268  0.000 0.000 6.50000
node 269  0.000 0.000 6.00000
node 270  0.000 0.000 5.50000
node 271  0.000 0.000 5.00000
node 272  0.000 0.000 4.50000
node 273  0.000 0.000 4.00000
node 274  0.000 0.000 3.50000
node 275  0.000 0.000 3.00000
node 276  0.000 0.000 2.50000
node 277  0.000 0.000 2.00000
node 278  0.000 0.000 1.50000
node 279  0.000 0.000 1.00000
node 280  0.000 0.000 0.50000
node 281  0.000 0.000 0.00000
puts "Finshed creating all p-y spring nodes..."

# define fixities for the soil end of springs
fix 2 1 1 1
fix 3 1 1 1
fix 4 1 1 1
fix 5 1 1 1
fix 6 1 1 1
fix 7 1 1 1
fix 8 1 1 1
fix 9 1 1 1
fix 10 1 1 1
fix 11 1 1 1
fix 12 1 1 1
fix 13 1 1 1
fix 14 1 1 1
fix 15 1 1 1
fix 16 1 1 1
fix 17 1 1 1
fix 18 1 1 1
fix 19 1 1 1
fix 20 1 1 1
fix 21 1 1 1
fix 22 1 1 1
fix 23 1 1 1
fix 24 1 1 1
fix 25 1 1 1
fix 26 1 1 1
fix 27 1 1 1
fix 28 1 1 1
fix 29 1 1 1
fix 30 1 1 1
fix 31 1 1 1
fix 32 1 1 1
fix 33 1 1 1
fix 34 1 1 1
fix 35 1 1 1
fix 36 1 1 1
fix 37 1 1 1
fix 38 1 1 1
fix 39 1 1 1
fix 40 1 1 1
fix 41 1 1 1
fix 42 1 1 1
fix 43 1 1 1
fix 44 1 1 1
fix 45 1 1 1
fix 46 1 1 1
fix 47 1 1 1
fix 48 1 1 1
fix 49 1 1 1
fix 50 1 1 1
fix 51 1 1 1
fix 52 1 1 1
fix 53 1 1 1
fix 54 1 1 1
fix 55 1 1 1
fix 56 1 1 1
fix 57 1 1 1
fix 58 1 1 1
fix 59 1 1 1
fix 60 1 1 1
fix 61 1 1 1
fix 62 1 1 1
fix 63 1 1 1
fix 64 1 1 1
fix 65 1 1 1
fix 66 1 1 1
fix 67 1 1 1
fix 68 1 1 1
fix 69 1 1 1
fix 70 1 1 1
fix 71 1 1 1
fix 72 1 1 1
fix 73 1 1 1
fix 74 1 1 1
fix 75 1 1 1
fix 76 1 1 1
fix 77 1 1 1
fix 78 1 1 1
fix 79 1 1 1
fix 80 1 1 1
fix 81 1 1 1

# define fixities for pile end of springs
fix 202 0 1 1
fix 203 0 1 1
fix 204 0 1 1
fix 205 0 1 1
fix 206 0 1 1
fix 207 0 1 1
fix 208 0 1 1
fix 209 0 1 1
fix 210 0 1 1
fix 211 0 1 1
fix 212 0 1 1
fix 213 0 1 1
fix 214 0 1 1
fix 215 0 1 1
fix 216 0 1 1
fix 217 0 1 1
fix 218 0 1 1
fix 219 0 1 1
fix 220 0 1 1
fix 221 0 1 1
fix 222 0 1 1
fix 223 0 1 1
fix 224 0 1 1
fix 225 0 1 1
fix 226 0 1 1
fix 227 0 1 1
fix 228 0 1 1
fix 229 0 1 1
fix 230 0 1 1
fix 231 0 1 1
fix 232 0 1 1
fix 233 0 1 1
fix 234 0 1 1
fix 235 0 1 1
fix 236 0 1 1
fix 237 0 1 1
fix 238 0 1 1
fix 239 0 1 1
fix 240 0 1 1
fix 241 0 1 1
fix 242 0 1 1
fix 243 0 1 1
fix 244 0 1 1
fix 245 0 1 1
fix 246 0 1 1
fix 247 0 1 1
fix 248 0 1 1
fix 249 0 1 1
fix 250 0 1 1
fix 251 0 1 1
fix 252 0 1 1
fix 253 0 1 1
fix 254 0 1 1
fix 255 0 1 1
fix 256 0 1 1
fix 257 0 1 1
fix 258 0 1 1
fix 259 0 1 1
fix 260 0 1 1
fix 261 0 1 1
fix 262 0 1 1
fix 263 0 1 1
fix 264 0 1 1
fix 265 0 1 1
fix 266 0 1 1
fix 267 0 1 1
fix 268 0 1 1
fix 269 0 1 1
fix 270 0 1 1
fix 271 0 1 1
fix 272 0 1 1
fix 273 0 1 1
fix 274 0 1 1
fix 275 0 1 1
fix 276 0 1 1
fix 277 0 1 1
fix 278 0 1 1
fix 279 0 1 1
fix 280 0 1 1
fix 281 0 1 1
puts "Finished creating all p-y spring boundary conditions..."

#-----------------------------------------------------------------------------------------
#  2. CREATE P-Y SPRING MATERIAL
#-----------------------------------------------------------------------------------------

uniaxialMaterial PySimple1 2 2     23.66983      0.00083 0
uniaxialMaterial PySimple1 3 2     54.93742      0.00136 0
uniaxialMaterial PySimple1 4 2     92.76801      0.00188 0
uniaxialMaterial PySimple1 5 2    136.28632      0.00239 0
uniaxialMaterial PySimple1 6 2    184.72733      0.00290 0
uniaxialMaterial PySimple1 7 2    237.38419      0.00340 0
uniaxialMaterial PySimple1 8 2    293.54190      0.00390 0
uniaxialMaterial PySimple1 9 2    352.37876      0.00438 0
uniaxialMaterial PySimple1 10 2    412.80561      0.00483 0
uniaxialMaterial PySimple1 11 2    473.19086      0.00526 0
uniaxialMaterial PySimple1 12 2    530.88085      0.00562 0
uniaxialMaterial PySimple1 13 2    581.35709      0.00590 0
uniaxialMaterial PySimple1 14 2    616.75420      0.00603 0
uniaxialMaterial PySimple1 15 2    623.25659      0.00594 0
uniaxialMaterial PySimple1 16 2    576.53590      0.00552 0
uniaxialMaterial PySimple1 17 2      9.30787      0.00082 0
uniaxialMaterial PySimple1 18 2     10.10868      0.00086 0
uniaxialMaterial PySimple1 19 2     10.92285      0.00090 0
uniaxialMaterial PySimple1 20 2     11.74934      0.00095 0
uniaxialMaterial PySimple1 21 2     12.58721      0.00099 0
uniaxialMaterial PySimple1 22 2    718.23881      0.00557 0
uniaxialMaterial PySimple1 23 2    984.69151      0.00738 0
uniaxialMaterial PySimple1 24 2   1200.82443      0.00880 0
uniaxialMaterial PySimple1 25 2   1380.75888      0.00990 0
uniaxialMaterial PySimple1 26 2   1534.82055      0.01078 0
uniaxialMaterial PySimple1 27 2   1670.51070      0.01151 0
uniaxialMaterial PySimple1 28 2   1793.24680      0.01212 0
uniaxialMaterial PySimple1 29 2   1906.92044      0.01266 0
uniaxialMaterial PySimple1 30 2   2014.31305      0.01314 0
uniaxialMaterial PySimple1 31 2   2117.40302      0.01358 0
uniaxialMaterial PySimple1 32 2   2217.59083      0.01399 0
uniaxialMaterial PySimple1 33 2   2315.86293      0.01438 0
uniaxialMaterial PySimple1 34 2   2412.91052      0.01476 0
uniaxialMaterial PySimple1 35 2   2509.21504      0.01512 0
uniaxialMaterial PySimple1 36 2   2605.10947      0.01547 0
uniaxialMaterial PySimple1 37 2   2700.82215      0.01581 0
uniaxialMaterial PySimple1 38 2   2796.50793      0.01615 0
uniaxialMaterial PySimple1 39 2   2892.27025      0.01648 0
uniaxialMaterial PySimple1 40 2   2988.17673      0.01681 0
uniaxialMaterial PySimple1 41 2   3084.27020      0.01713 0
uniaxialMaterial PySimple1 42 2   3180.57640      0.01745 0
uniaxialMaterial PySimple1 43 2   3277.10936      0.01776 0
uniaxialMaterial PySimple1 44 2   3373.87523      0.01808 0
uniaxialMaterial PySimple1 45 2   3470.87482      0.01838 0
uniaxialMaterial PySimple1 46 2   3568.10549      0.01869 0
uniaxialMaterial PySimple1 47 2   3665.56234      0.01899 0
uniaxialMaterial PySimple1 48 2   3763.23912      0.01928 0
uniaxialMaterial PySimple1 49 2   3861.12878      0.01958 0
uniaxialMaterial PySimple1 50 2   3959.22389      0.01987 0
uniaxialMaterial PySimple1 51 2   4057.51687      0.02016 0
uniaxialMaterial PySimple1 52 2   4156.00020      0.02044 0
uniaxialMaterial PySimple1 53 2   4254.66650      0.02073 0
uniaxialMaterial PySimple1 54 2   4353.50862      0.02101 0
uniaxialMaterial PySimple1 55 2   4452.51967      0.02129 0
uniaxialMaterial PySimple1 56 2   4551.69304      0.02156 0
uniaxialMaterial PySimple1 57 2   4651.02242      0.02183 0
uniaxialMaterial PySimple1 58 2   4750.50177      0.02210 0
uniaxialMaterial PySimple1 59 2   4850.12535      0.02237 0
uniaxialMaterial PySimple1 60 2   4949.88767      0.02264 0
uniaxialMaterial PySimple1 61 2   5049.78354      0.02290 0
uniaxialMaterial PySimple1 62 2   5149.80798      0.02316 0
uniaxialMaterial PySimple1 63 2   5249.95628      0.02342 0
uniaxialMaterial PySimple1 64 2   5350.22393      0.02368 0
uniaxialMaterial PySimple1 65 2   5450.60664      0.02394 0
uniaxialMaterial PySimple1 66 2   5551.10034      0.02419 0
uniaxialMaterial PySimple1 67 2   5651.70112      0.02444 0
uniaxialMaterial PySimple1 68 2   5752.40527      0.02469 0
uniaxialMaterial PySimple1 69 2   5853.20924      0.02494 0
uniaxialMaterial PySimple1 70 2   5954.10964      0.02518 0
uniaxialMaterial PySimple1 71 2   6055.10324      0.02542 0
uniaxialMaterial PySimple1 72 2   6156.18694      0.02567 0
uniaxialMaterial PySimple1 73 2   6257.35780      0.02591 0
uniaxialMaterial PySimple1 74 2   6358.61297      0.02614 0
uniaxialMaterial PySimple1 75 2   6459.94975      0.02638 0
uniaxialMaterial PySimple1 76 2   6561.36555      0.02662 0
uniaxialMaterial PySimple1 77 2   6662.85788      0.02685 0
uniaxialMaterial PySimple1 78 2   6764.42436      0.02708 0
uniaxialMaterial PySimple1 79 2   6866.06270      0.02731 0
uniaxialMaterial PySimple1 80 2   6967.77072      0.02754 0
uniaxialMaterial PySimple1 81 2   7069.54631      0.02777 0
puts "Finished creating all p-y spring materials..."

#-----------------------------------------------------------------------------------------
#  3. CREATE P-Y SPRING ELEMENTS
#-----------------------------------------------------------------------------------------

element zeroLength 2 2 202 -mat 2 -dir 1
element zeroLength 3 3 203 -mat 3 -dir 1
element zeroLength 4 4 204 -mat 4 -dir 1
element zeroLength 5 5 205 -mat 5 -dir 1
element zeroLength 6 6 206 -mat 6 -dir 1
element zeroLength 7 7 207 -mat 7 -dir 1
element zeroLength 8 8 208 -mat 8 -dir 1
element zeroLength 9 9 209 -mat 9 -dir 1
element zeroLength 10 10 210 -mat 10 -dir 1
element zeroLength 11 11 211 -mat 11 -dir 1
element zeroLength 12 12 212 -mat 12 -dir 1
element zeroLength 13 13 213 -mat 13 -dir 1
element zeroLength 14 14 214 -mat 14 -dir 1
element zeroLength 15 15 215 -mat 15 -dir 1
element zeroLength 16 16 216 -mat 16 -dir 1
element zeroLength 17 17 217 -mat 17 -dir 1
element zeroLength 18 18 218 -mat 18 -dir 1
element zeroLength 19 19 219 -mat 19 -dir 1
element zeroLength 20 20 220 -mat 20 -dir 1
element zeroLength 21 21 221 -mat 21 -dir 1
element zeroLength 22 22 222 -mat 22 -dir 1
element zeroLength 23 23 223 -mat 23 -dir 1
element zeroLength 24 24 224 -mat 24 -dir 1
element zeroLength 25 25 225 -mat 25 -dir 1
element zeroLength 26 26 226 -mat 26 -dir 1
element zeroLength 27 27 227 -mat 27 -dir 1
element zeroLength 28 28 228 -mat 28 -dir 1
element zeroLength 29 29 229 -mat 29 -dir 1
element zeroLength 30 30 230 -mat 30 -dir 1
element zeroLength 31 31 231 -mat 31 -dir 1
element zeroLength 32 32 232 -mat 32 -dir 1
element zeroLength 33 33 233 -mat 33 -dir 1
element zeroLength 34 34 234 -mat 34 -dir 1
element zeroLength 35 35 235 -mat 35 -dir 1
element zeroLength 36 36 236 -mat 36 -dir 1
element zeroLength 37 37 237 -mat 37 -dir 1
element zeroLength 38 38 238 -mat 38 -dir 1
element zeroLength 39 39 239 -mat 39 -dir 1
element zeroLength 40 40 240 -mat 40 -dir 1
element zeroLength 41 41 241 -mat 41 -dir 1
element zeroLength 42 42 242 -mat 42 -dir 1
element zeroLength 43 43 243 -mat 43 -dir 1
element zeroLength 44 44 244 -mat 44 -dir 1
element zeroLength 45 45 245 -mat 45 -dir 1
element zeroLength 46 46 246 -mat 46 -dir 1
element zeroLength 47 47 247 -mat 47 -dir 1
element zeroLength 48 48 248 -mat 48 -dir 1
element zeroLength 49 49 249 -mat 49 -dir 1
element zeroLength 50 50 250 -mat 50 -dir 1
element zeroLength 51 51 251 -mat 51 -dir 1
element zeroLength 52 52 252 -mat 52 -dir 1
element zeroLength 53 53 253 -mat 53 -dir 1
element zeroLength 54 54 254 -mat 54 -dir 1
element zeroLength 55 55 255 -mat 55 -dir 1
element zeroLength 56 56 256 -mat 56 -dir 1
element zeroLength 57 57 257 -mat 57 -dir 1
element zeroLength 58 58 258 -mat 58 -dir 1
element zeroLength 59 59 259 -mat 59 -dir 1
element zeroLength 60 60 260 -mat 60 -dir 1
element zeroLength 61 61 261 -mat 61 -dir 1
element zeroLength 62 62 262 -mat 62 -dir 1
element zeroLength 63 63 263 -mat 63 -dir 1
element zeroLength 64 64 264 -mat 64 -dir 1
element zeroLength 65 65 265 -mat 65 -dir 1
element zeroLength 66 66 266 -mat 66 -dir 1
element zeroLength 67 67 267 -mat 67 -dir 1
element zeroLength 68 68 268 -mat 68 -dir 1
element zeroLength 69 69 269 -mat 69 -dir 1
element zeroLength 70 70 270 -mat 70 -dir 1
element zeroLength 71 71 271 -mat 71 -dir 1
element zeroLength 72 72 272 -mat 72 -dir 1
element zeroLength 73 73 273 -mat 73 -dir 1
element zeroLength 74 74 274 -mat 74 -dir 1
element zeroLength 75 75 275 -mat 75 -dir 1
element zeroLength 76 76 276 -mat 76 -dir 1
element zeroLength 77 77 277 -mat 77 -dir 1
element zeroLength 78 78 278 -mat 78 -dir 1
element zeroLength 79 79 279 -mat 79 -dir 1
element zeroLength 80 80 280 -mat 80 -dir 1
element zeroLength 81 81 281 -mat 81 -dir 1
puts "Finished creating all p-y spring elements..."

#-----------------------------------------------------------------------------------------
#  4. CREATE PILE NODES AND BOUNDARY CONDITIONS
#-----------------------------------------------------------------------------------------
model BasicBuilder -ndm 3 -ndf 6

# create nodes
node 501  0.000 0.000 40.00000
node 502  0.000 0.000 39.50000
node 503  0.000 0.000 39.00000
node 504  0.000 0.000 38.50000
node 505  0.000 0.000 38.00000
node 506  0.000 0.000 37.50000
node 507  0.000 0.000 37.00000
node 508  0.000 0.000 36.50000
node 509  0.000 0.000 36.00000
node 510  0.000 0.000 35.50000
node 511  0.000 0.000 35.00000
node 512  0.000 0.000 34.50000
node 513  0.000 0.000 34.00000
node 514  0.000 0.000 33.50000
node 515  0.000 0.000 33.00000
node 516  0.000 0.000 32.50000
node 517  0.000 0.000 32.00000
node 518  0.000 0.000 31.50000
node 519  0.000 0.000 31.00000
node 520  0.000 0.000 30.50000
node 521  0.000 0.000 30.00000
node 522  0.000 0.000 29.50000
node 523  0.000 0.000 29.00000
node 524  0.000 0.000 28.50000
node 525  0.000 0.000 28.00000
node 526  0.000 0.000 27.50000
node 527  0.000 0.000 27.00000
node 528  0.000 0.000 26.50000
node 529  0.000 0.000 26.00000
node 530  0.000 0.000 25.50000
node 531  0.000 0.000 25.00000
node 532  0.000 0.000 24.50000
node 533  0.000 0.000 24.00000
node 534  0.000 0.000 23.50000
node 535  0.000 0.000 23.00000
node 536  0.000 0.000 22.50000
node 537  0.000 0.000 22.00000
node 538  0.000 0.000 21.50000
node 539  0.000 0.000 21.00000
node 540  0.000 0.000 20.50000
node 541  0.000 0.000 20.00000
node 542  0.000 0.000 19.50000
node 543  0.000 0.000 19.00000
node 544  0.000 0.000 18.50000
node 545  0.000 0.000 18.00000
node 546  0.000 0.000 17.50000
node 547  0.000 0.000 17.00000
node 548  0.000 0.000 16.50000
node 549  0.000 0.000 16.00000
node 550  0.000 0.000 15.50000
node 551  0.000 0.000 15.00000
node 552  0.000 0.000 14.50000
node 553  0.000 0.000 14.00000
node 554  0.000 0.000 13.50000
node 555  0.000 0.000 13.00000
node 556  0.000 0.000 12.50000
node 557  0.000 0.000 12.00000
node 558  0.000 0.000 11.50000
node 559  0.000 0.000 11.00000
node 560  0.000 0.000 10.50000
node 561  0.000 0.000 10.00000
node 562  0.000 0.000 9.50000
node 563  0.000 0.000 9.00000
node 564  0.000 0.000 8.50000
node 565  0.000 0.000 8.00000
node 566  0.000 0.000 7.50000
node 567  0.000 0.000 7.00000
node 568  0.000 0.000 6.50000
node 569  0.000 0.000 6.00000
node 570  0.000 0.000 5.50000
node 571  0.000 0.000 5.00000
node 572  0.000 0.000 4.50000
node 573  0.000 0.000 4.00000
node 574  0.000 0.000 3.50000
node 575  0.000 0.000 3.00000
node 576  0.000 0.000 2.50000
node 577  0.000 0.000 2.00000
node 578  0.000 0.000 1.50000
node 579  0.000 0.000 1.00000
node 580  0.000 0.000 0.50000
node 581  0.000 0.000 0.00000

set nodesInfo6 [open $outDir/NodesInfo6.dat w]
puts $nodesInfo6 "501 0 0 40.00000"
puts $nodesInfo6 "502 0 0 39.50000"
puts $nodesInfo6 "503 0 0 39.00000"
puts $nodesInfo6 "504 0 0 38.50000"
puts $nodesInfo6 "505 0 0 38.00000"
puts $nodesInfo6 "506 0 0 37.50000"
puts $nodesInfo6 "507 0 0 37.00000"
puts $nodesInfo6 "508 0 0 36.50000"
puts $nodesInfo6 "509 0 0 36.00000"
puts $nodesInfo6 "510 0 0 35.50000"
puts $nodesInfo6 "511 0 0 35.00000"
puts $nodesInfo6 "512 0 0 34.50000"
puts $nodesInfo6 "513 0 0 34.00000"
puts $nodesInfo6 "514 0 0 33.50000"
puts $nodesInfo6 "515 0 0 33.00000"
puts $nodesInfo6 "516 0 0 32.50000"
puts $nodesInfo6 "517 0 0 32.00000"
puts $nodesInfo6 "518 0 0 31.50000"
puts $nodesInfo6 "519 0 0 31.00000"
puts $nodesInfo6 "520 0 0 30.50000"
puts $nodesInfo6 "521 0 0 30.00000"
puts $nodesInfo6 "522 0 0 29.50000"
puts $nodesInfo6 "523 0 0 29.00000"
puts $nodesInfo6 "524 0 0 28.50000"
puts $nodesInfo6 "525 0 0 28.00000"
puts $nodesInfo6 "526 0 0 27.50000"
puts $nodesInfo6 "527 0 0 27.00000"
puts $nodesInfo6 "528 0 0 26.50000"
puts $nodesInfo6 "529 0 0 26.00000"
puts $nodesInfo6 "530 0 0 25.50000"
puts $nodesInfo6 "531 0 0 25.00000"
puts $nodesInfo6 "532 0 0 24.50000"
puts $nodesInfo6 "533 0 0 24.00000"
puts $nodesInfo6 "534 0 0 23.50000"
puts $nodesInfo6 "535 0 0 23.00000"
puts $nodesInfo6 "536 0 0 22.50000"
puts $nodesInfo6 "537 0 0 22.00000"
puts $nodesInfo6 "538 0 0 21.50000"
puts $nodesInfo6 "539 0 0 21.00000"
puts $nodesInfo6 "540 0 0 20.50000"
puts $nodesInfo6 "541 0 0 20.00000"
puts $nodesInfo6 "542 0 0 19.50000"
puts $nodesInfo6 "543 0 0 19.00000"
puts $nodesInfo6 "544 0 0 18.50000"
puts $nodesInfo6 "545 0 0 18.00000"
puts $nodesInfo6 "546 0 0 17.50000"
puts $nodesInfo6 "547 0 0 17.00000"
puts $nodesInfo6 "548 0 0 16.50000"
puts $nodesInfo6 "549 0 0 16.00000"
puts $nodesInfo6 "550 0 0 15.50000"
puts $nodesInfo6 "551 0 0 15.00000"
puts $nodesInfo6 "552 0 0 14.50000"
puts $nodesInfo6 "553 0 0 14.00000"
puts $nodesInfo6 "554 0 0 13.50000"
puts $nodesInfo6 "555 0 0 13.00000"
puts $nodesInfo6 "556 0 0 12.50000"
puts $nodesInfo6 "557 0 0 12.00000"
puts $nodesInfo6 "558 0 0 11.50000"
puts $nodesInfo6 "559 0 0 11.00000"
puts $nodesInfo6 "560 0 0 10.50000"
puts $nodesInfo6 "561 0 0 10.00000"
puts $nodesInfo6 "562 0 0 9.50000"
puts $nodesInfo6 "563 0 0 9.00000"
puts $nodesInfo6 "564 0 0 8.50000"
puts $nodesInfo6 "565 0 0 8.00000"
puts $nodesInfo6 "566 0 0 7.50000"
puts $nodesInfo6 "567 0 0 7.00000"
puts $nodesInfo6 "568 0 0 6.50000"
puts $nodesInfo6 "569 0 0 6.00000"
puts $nodesInfo6 "570 0 0 5.50000"
puts $nodesInfo6 "571 0 0 5.00000"
puts $nodesInfo6 "572 0 0 4.50000"
puts $nodesInfo6 "573 0 0 4.00000"
puts $nodesInfo6 "574 0 0 3.50000"
puts $nodesInfo6 "575 0 0 3.00000"
puts $nodesInfo6 "576 0 0 2.50000"
puts $nodesInfo6 "577 0 0 2.00000"
puts $nodesInfo6 "578 0 0 1.50000"
puts $nodesInfo6 "579 0 0 1.00000"
puts $nodesInfo6 "580 0 0 0.50000"
puts $nodesInfo6 "581 0 0 0.00000"
close $nodesInfo6
puts "Finished creating all pile nodes..."

# create geometric transformation for vertical pile
set transTag 1
geomTransf Linear $transTag  0 -1 0

# create the fixities for the pile nodes
fix 501  0 1 0 1 0 1
fix 502  0 1 0 1 0 1
fix 503  0 1 0 1 0 1
fix 504  0 1 0 1 0 1
fix 505  0 1 0 1 0 1
fix 506  0 1 0 1 0 1
fix 507  0 1 0 1 0 1
fix 508  0 1 0 1 0 1
fix 509  0 1 0 1 0 1
fix 510  0 1 0 1 0 1
fix 511  0 1 0 1 0 1
fix 512  0 1 0 1 0 1
fix 513  0 1 0 1 0 1
fix 514  0 1 0 1 0 1
fix 515  0 1 0 1 0 1
fix 516  0 1 0 1 0 1
fix 517  0 1 0 1 0 1
fix 518  0 1 0 1 0 1
fix 519  0 1 0 1 0 1
fix 520  0 1 0 1 0 1
fix 521  0 1 0 1 0 1
fix 522  0 1 0 1 0 1
fix 523  0 1 0 1 0 1
fix 524  0 1 0 1 0 1
fix 525  0 1 0 1 0 1
fix 526  0 1 0 1 0 1
fix 527  0 1 0 1 0 1
fix 528  0 1 0 1 0 1
fix 529  0 1 0 1 0 1
fix 530  0 1 0 1 0 1
fix 531  0 1 0 1 0 1
fix 532  0 1 0 1 0 1
fix 533  0 1 0 1 0 1
fix 534  0 1 0 1 0 1
fix 535  0 1 0 1 0 1
fix 536  0 1 0 1 0 1
fix 537  0 1 0 1 0 1
fix 538  0 1 0 1 0 1
fix 539  0 1 0 1 0 1
fix 540  0 1 0 1 0 1
fix 541  0 1 0 1 0 1
fix 542  0 1 0 1 0 1
fix 543  0 1 0 1 0 1
fix 544  0 1 0 1 0 1
fix 545  0 1 0 1 0 1
fix 546  0 1 0 1 0 1
fix 547  0 1 0 1 0 1
fix 548  0 1 0 1 0 1
fix 549  0 1 0 1 0 1
fix 550  0 1 0 1 0 1
fix 551  0 1 0 1 0 1
fix 552  0 1 0 1 0 1
fix 553  0 1 0 1 0 1
fix 554  0 1 0 1 0 1
fix 555  0 1 0 1 0 1
fix 556  0 1 0 1 0 1
fix 557  0 1 0 1 0 1
fix 558  0 1 0 1 0 1
fix 559  0 1 0 1 0 1
fix 560  0 1 0 1 0 1
fix 561  0 1 0 1 0 1
fix 562  0 1 0 1 0 1
fix 563  0 1 0 1 0 1
fix 564  0 1 0 1 0 1
fix 565  0 1 0 1 0 1
fix 566  0 1 0 1 0 1
fix 567  0 1 0 1 0 1
fix 568  0 1 0 1 0 1
fix 569  0 1 0 1 0 1
fix 570  0 1 0 1 0 1
fix 571  0 1 0 1 0 1
fix 572  0 1 0 1 0 1
fix 573  0 1 0 1 0 1
fix 574  0 1 0 1 0 1
fix 575  0 1 0 1 0 1
fix 576  0 1 0 1 0 1
fix 577  0 1 0 1 0 1
fix 578  0 1 0 1 0 1
fix 579  0 1 0 1 0 1
fix 580  0 1 0 1 0 1
fix 581  0 1 1 1 0 1
puts "Finished creating all pile boundary conditions..."

# define equalDOF for pile and spring nodes
equalDOF 502 202 1
equalDOF 503 203 1
equalDOF 504 204 1
equalDOF 505 205 1
equalDOF 506 206 1
equalDOF 507 207 1
equalDOF 508 208 1
equalDOF 509 209 1
equalDOF 510 210 1
equalDOF 511 211 1
equalDOF 512 212 1
equalDOF 513 213 1
equalDOF 514 214 1
equalDOF 515 215 1
equalDOF 516 216 1
equalDOF 517 217 1
equalDOF 518 218 1
equalDOF 519 219 1
equalDOF 520 220 1
equalDOF 521 221 1
equalDOF 522 222 1
equalDOF 523 223 1
equalDOF 524 224 1
equalDOF 525 225 1
equalDOF 526 226 1
equalDOF 527 227 1
equalDOF 528 228 1
equalDOF 529 229 1
equalDOF 530 230 1
equalDOF 531 231 1
equalDOF 532 232 1
equalDOF 533 233 1
equalDOF 534 234 1
equalDOF 535 235 1
equalDOF 536 236 1
equalDOF 537 237 1
equalDOF 538 238 1
equalDOF 539 239 1
equalDOF 540 240 1
equalDOF 541 241 1
equalDOF 542 242 1
equalDOF 543 243 1
equalDOF 544 244 1
equalDOF 545 245 1
equalDOF 546 246 1
equalDOF 547 247 1
equalDOF 548 248 1
equalDOF 549 249 1
equalDOF 550 250 1
equalDOF 551 251 1
equalDOF 552 252 1
equalDOF 553 253 1
equalDOF 554 254 1
equalDOF 555 255 1
equalDOF 556 256 1
equalDOF 557 257 1
equalDOF 558 258 1
equalDOF 559 259 1
equalDOF 560 260 1
equalDOF 561 261 1
equalDOF 562 262 1
equalDOF 563 263 1
equalDOF 564 264 1
equalDOF 565 265 1
equalDOF 566 266 1
equalDOF 567 267 1
equalDOF 568 268 1
equalDOF 569 269 1
equalDOF 570 270 1
equalDOF 571 271 1
equalDOF 572 272 1
equalDOF 573 273 1
equalDOF 574 274 1
equalDOF 575 275 1
equalDOF 576 276 1
equalDOF 577 277 1
equalDOF 578 278 1
equalDOF 579 279 1
equalDOF 580 280 1
equalDOF 581 281 1
puts "Finished creating all pile-spring equalDOF..."

#-----------------------------------------------------------------------------------------
#  5. CREATE ELASTIC PILE SECTION
#-----------------------------------------------------------------------------------------

# pile properties
set pi       [expr 4.0*atan(1.0)]
set radius   0.30000
set area     [expr $pi*$radius*$radius]
set I        [expr $pi*$radius*$radius*$radius*$radius/4.0]
set J        [expr 2.0*$I]
set E        30000000.00000
set nu       0.30000
set G        [expr $E/(2.0*(1.0+$nu))]
set nIntPts  3
set secTag   1
set secTag3D 3

# beam section, no torsion
section Elastic  $secTag $E $area $I $I $G $J
# create torsional section
uniaxialMaterial Elastic 100 1.e10
# aggregate the beam section
section Aggregator $secTag3D  100 T -section $secTag
puts "Finished creating all beam materials..."

#-----------------------------------------------------------------------------------------
#  6. CREATE PILE ELEMENTS
#-----------------------------------------------------------------------------------------

element dispBeamColumn 501 501 502 $nIntPts $secTag3D $transTag
element dispBeamColumn 502 502 503 $nIntPts $secTag3D $transTag
element dispBeamColumn 503 503 504 $nIntPts $secTag3D $transTag
element dispBeamColumn 504 504 505 $nIntPts $secTag3D $transTag
element dispBeamColumn 505 505 506 $nIntPts $secTag3D $transTag
element dispBeamColumn 506 506 507 $nIntPts $secTag3D $transTag
element dispBeamColumn 507 507 508 $nIntPts $secTag3D $transTag
element dispBeamColumn 508 508 509 $nIntPts $secTag3D $transTag
element dispBeamColumn 509 509 510 $nIntPts $secTag3D $transTag
element dispBeamColumn 510 510 511 $nIntPts $secTag3D $transTag
element dispBeamColumn 511 511 512 $nIntPts $secTag3D $transTag
element dispBeamColumn 512 512 513 $nIntPts $secTag3D $transTag
element dispBeamColumn 513 513 514 $nIntPts $secTag3D $transTag
element dispBeamColumn 514 514 515 $nIntPts $secTag3D $transTag
element dispBeamColumn 515 515 516 $nIntPts $secTag3D $transTag
element dispBeamColumn 516 516 517 $nIntPts $secTag3D $transTag
element dispBeamColumn 517 517 518 $nIntPts $secTag3D $transTag
element dispBeamColumn 518 518 519 $nIntPts $secTag3D $transTag
element dispBeamColumn 519 519 520 $nIntPts $secTag3D $transTag
element dispBeamColumn 520 520 521 $nIntPts $secTag3D $transTag
element dispBeamColumn 521 521 522 $nIntPts $secTag3D $transTag
element dispBeamColumn 522 522 523 $nIntPts $secTag3D $transTag
element dispBeamColumn 523 523 524 $nIntPts $secTag3D $transTag
element dispBeamColumn 524 524 525 $nIntPts $secTag3D $transTag
element dispBeamColumn 525 525 526 $nIntPts $secTag3D $transTag
element dispBeamColumn 526 526 527 $nIntPts $secTag3D $transTag
element dispBeamColumn 527 527 528 $nIntPts $secTag3D $transTag
element dispBeamColumn 528 528 529 $nIntPts $secTag3D $transTag
element dispBeamColumn 529 529 530 $nIntPts $secTag3D $transTag
element dispBeamColumn 530 530 531 $nIntPts $secTag3D $transTag
element dispBeamColumn 531 531 532 $nIntPts $secTag3D $transTag
element dispBeamColumn 532 532 533 $nIntPts $secTag3D $transTag
element dispBeamColumn 533 533 534 $nIntPts $secTag3D $transTag
element dispBeamColumn 534 534 535 $nIntPts $secTag3D $transTag
element dispBeamColumn 535 535 536 $nIntPts $secTag3D $transTag
element dispBeamColumn 536 536 537 $nIntPts $secTag3D $transTag
element dispBeamColumn 537 537 538 $nIntPts $secTag3D $transTag
element dispBeamColumn 538 538 539 $nIntPts $secTag3D $transTag
element dispBeamColumn 539 539 540 $nIntPts $secTag3D $transTag
element dispBeamColumn 540 540 541 $nIntPts $secTag3D $transTag
element dispBeamColumn 541 541 542 $nIntPts $secTag3D $transTag
element dispBeamColumn 542 542 543 $nIntPts $secTag3D $transTag
element dispBeamColumn 543 543 544 $nIntPts $secTag3D $transTag
element dispBeamColumn 544 544 545 $nIntPts $secTag3D $transTag
element dispBeamColumn 545 545 546 $nIntPts $secTag3D $transTag
element dispBeamColumn 546 546 547 $nIntPts $secTag3D $transTag
element dispBeamColumn 547 547 548 $nIntPts $secTag3D $transTag
element dispBeamColumn 548 548 549 $nIntPts $secTag3D $transTag
element dispBeamColumn 549 549 550 $nIntPts $secTag3D $transTag
element dispBeamColumn 550 550 551 $nIntPts $secTag3D $transTag
element dispBeamColumn 551 551 552 $nIntPts $secTag3D $transTag
element dispBeamColumn 552 552 553 $nIntPts $secTag3D $transTag
element dispBeamColumn 553 553 554 $nIntPts $secTag3D $transTag
element dispBeamColumn 554 554 555 $nIntPts $secTag3D $transTag
element dispBeamColumn 555 555 556 $nIntPts $secTag3D $transTag
element dispBeamColumn 556 556 557 $nIntPts $secTag3D $transTag
element dispBeamColumn 557 557 558 $nIntPts $secTag3D $transTag
element dispBeamColumn 558 558 559 $nIntPts $secTag3D $transTag
element dispBeamColumn 559 559 560 $nIntPts $secTag3D $transTag
element dispBeamColumn 560 560 561 $nIntPts $secTag3D $transTag
element dispBeamColumn 561 561 562 $nIntPts $secTag3D $transTag
element dispBeamColumn 562 562 563 $nIntPts $secTag3D $transTag
element dispBeamColumn 563 563 564 $nIntPts $secTag3D $transTag
element dispBeamColumn 564 564 565 $nIntPts $secTag3D $transTag
element dispBeamColumn 565 565 566 $nIntPts $secTag3D $transTag
element dispBeamColumn 566 566 567 $nIntPts $secTag3D $transTag
element dispBeamColumn 567 567 568 $nIntPts $secTag3D $transTag
element dispBeamColumn 568 568 569 $nIntPts $secTag3D $transTag
element dispBeamColumn 569 569 570 $nIntPts $secTag3D $transTag
element dispBeamColumn 570 570 571 $nIntPts $secTag3D $transTag
element dispBeamColumn 571 571 572 $nIntPts $secTag3D $transTag
element dispBeamColumn 572 572 573 $nIntPts $secTag3D $transTag
element dispBeamColumn 573 573 574 $nIntPts $secTag3D $transTag
element dispBeamColumn 574 574 575 $nIntPts $secTag3D $transTag
element dispBeamColumn 575 575 576 $nIntPts $secTag3D $transTag
element dispBeamColumn 576 576 577 $nIntPts $secTag3D $transTag
element dispBeamColumn 577 577 578 $nIntPts $secTag3D $transTag
element dispBeamColumn 578 578 579 $nIntPts $secTag3D $transTag
element dispBeamColumn 579 579 580 $nIntPts $secTag3D $transTag
element dispBeamColumn 580 580 581 $nIntPts $secTag3D $transTag

set pileInfo [open $outDir/PileInfo.dat w]
puts $pileInfo "501 501 502"
puts $pileInfo "502 502 503"
puts $pileInfo "503 503 504"
puts $pileInfo "504 504 505"
puts $pileInfo "505 505 506"
puts $pileInfo "506 506 507"
puts $pileInfo "507 507 508"
puts $pileInfo "508 508 509"
puts $pileInfo "509 509 510"
puts $pileInfo "510 510 511"
puts $pileInfo "511 511 512"
puts $pileInfo "512 512 513"
puts $pileInfo "513 513 514"
puts $pileInfo "514 514 515"
puts $pileInfo "515 515 516"
puts $pileInfo "516 516 517"
puts $pileInfo "517 517 518"
puts $pileInfo "518 518 519"
puts $pileInfo "519 519 520"
puts $pileInfo "520 520 521"
puts $pileInfo "521 521 522"
puts $pileInfo "522 522 523"
puts $pileInfo "523 523 524"
puts $pileInfo "524 524 525"
puts $pileInfo "525 525 526"
puts $pileInfo "526 526 527"
puts $pileInfo "527 527 528"
puts $pileInfo "528 528 529"
puts $pileInfo "529 529 530"
puts $pileInfo "530 530 531"
puts $pileInfo "531 531 532"
puts $pileInfo "532 532 533"
puts $pileInfo "533 533 534"
puts $pileInfo "534 534 535"
puts $pileInfo "535 535 536"
puts $pileInfo "536 536 537"
puts $pileInfo "537 537 538"
puts $pileInfo "538 538 539"
puts $pileInfo "539 539 540"
puts $pileInfo "540 540 541"
puts $pileInfo "541 541 542"
puts $pileInfo "542 542 543"
puts $pileInfo "543 543 544"
puts $pileInfo "544 544 545"
puts $pileInfo "545 545 546"
puts $pileInfo "546 546 547"
puts $pileInfo "547 547 548"
puts $pileInfo "548 548 549"
puts $pileInfo "549 549 550"
puts $pileInfo "550 550 551"
puts $pileInfo "551 551 552"
puts $pileInfo "552 552 553"
puts $pileInfo "553 553 554"
puts $pileInfo "554 554 555"
puts $pileInfo "555 555 556"
puts $pileInfo "556 556 557"
puts $pileInfo "557 557 558"
puts $pileInfo "558 558 559"
puts $pileInfo "559 559 560"
puts $pileInfo "560 560 561"
puts $pileInfo "561 561 562"
puts $pileInfo "562 562 563"
puts $pileInfo "563 563 564"
puts $pileInfo "564 564 565"
puts $pileInfo "565 565 566"
puts $pileInfo "566 566 567"
puts $pileInfo "567 567 568"
puts $pileInfo "568 568 569"
puts $pileInfo "569 569 570"
puts $pileInfo "570 570 571"
puts $pileInfo "571 571 572"
puts $pileInfo "572 572 573"
puts $pileInfo "573 573 574"
puts $pileInfo "574 574 575"
puts $pileInfo "575 575 576"
puts $pileInfo "576 576 577"
puts $pileInfo "577 577 578"
puts $pileInfo "578 578 579"
puts $pileInfo "579 579 580"
puts $pileInfo "580 580 581"
close $pileInfo
puts "Finished creating all beam elements..."

#-----------------------------------------------------------------------------------------
#  7. CREATE RECORDERS
#-----------------------------------------------------------------------------------------

# create list with pile node info
set nodeList6 {}
set channel [open $outDir/"NodesInfo6.dat" r]
set ctr 0;
foreach line [split [read -nonewline $channel] \n] {
  set ctr0 [expr $ctr+1];
  set lineData($ctr) $line
  set nodeNumber [lindex $lineData($ctr) 0]
  lappend nodeList6 $nodeNumber
}
close $channel

# create pile element list
set BeamElementList {}
set channel [open $outDir/"PileInfo.dat" r]
set ctr 0;
foreach line [split [read -nonewline $channel] \n] {
  set ctr0 [expr $ctr+1];
  set lineData($ctr) $line
  set elementNumber [lindex $lineData($ctr) 0]
  lappend BeamElementList $elementNumber
}
close $channel

# designate recorder timestep
set dt 0.5

recorder Node -file pyDisplace.out -time -nodeRange 202 281 -dof 1 -dT $dt disp
recorder Node -file pyForces.out   -time -nodeRange 2   81 -dof 1 -dT $dt reaction

eval "recorder Node -file Displacements.out -time -node $nodeList6 -dof 1 2 3 -dT $dt disp"
eval "recorder Node -file Reaction.out      -time -node $nodeList6 -dof 1 2 3 -dT $dt reaction "
eval "recorder Element -file globalForces.out    -time -ele $BeamElementList -dT $dt globalForce"

#-----------------------------------------------------------------------------------------
#  8. CREATE APPLIED DISPLACEMENT PROFILE
#-----------------------------------------------------------------------------------------

pattern Plain 10 {Series -time {0 10 10000} -values {0 1 1} -factor 1} {
    sp 2 1  0.300
    sp 3 1  0.300
    sp 4 1  0.300
    sp 5 1  0.300
    sp 6 1  0.300
    sp 7 1  0.300
    sp 8 1  0.300
    sp 9 1  0.300
    sp 10 1  0.300
    sp 11 1  0.300
    sp 12 1  0.300
    sp 13 1  0.300
    sp 14 1  0.300
    sp 15 1  0.300
    sp 16 1  0.300
    sp 17 1  0.225
    sp 18 1  0.150
    sp 19 1  0.075
    sp 20 1  0.000
}
puts "Finished adding all sp constraints for push-over analysis..."

#-----------------------------------------------------------------------------------------
#  8. CREATE AND EXECUTE THE ANALYSIS
#-----------------------------------------------------------------------------------------

integrator  LoadControl 0.05
numberer    RCM
system      BandGeneral
constraints Transformation
test        NormDispIncr 1e-5 20 1
algorithm   Newton
analysis    Static

set startT [clock seconds]

analyze     205

set endT [clock seconds]

puts "Load Application finished..."
puts "Analysis execution time: [expr $endT-$startT] seconds"

wipe
