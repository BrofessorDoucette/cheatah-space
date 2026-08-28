// Unit tests for space.irbem's IRBEM-named adapter layer — api.hpp.
//
// The layer has no physics of its own, so these tests ask three questions and nothing else:
//
//   1. **Does the name mean what IRBEM's name means?** Answered differentially, against golden
//      values harvested from the reference library as a BLACK BOX (dlopen + the documented C entry
//      points on /tmp/irbem-builds/libirbem-O2.so; see tools/oracle/). The literals below carry the
//      exact call that produced them. Nothing here links IRBEM, so the QA gate can run this file on
//      a machine that has never seen the reference.
//   2. **Does the adapter compose its legs correctly?** Answered structurally: every frame pair is
//      round-tripped, and the batch lanes are compared to the scalar lane with `==` rather than a
//      tolerance, because they are required to compute the identical expression rather than an
//      equivalent one.
//   3. **Does an unrepresentable request get REPORTED?** IRBEM prints `sysaxesOUT out of range !`
//      on stdout and returns baddata; this layer returns Status::DomainError and writes nothing.
//
// On the epoch sweep: the transforms here are parameterized by EPOCH, not by geomagnetic activity —
// no routine in api.hpp reads a maginput vector, because no external field model exists in this
// module yet. So the analogue of "sample storms, not just quiet times" is to sample the epoch
// axis hard: the goldens span 1965-2029, both solstices, both equinoxes, midnight, noon and the
// last second of a day, which is where the sidereal-time and dipole-drift terms move most.
#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <iterator>
#include <limits>
#include <numbers>
#include <span>
#include <vector>

#include "alloc_counter.hpp"
#include "space/irbem/api.hpp"

namespace ib = cheatah::space::irbem;
namespace api = ib::api;
namespace fx = cheatah::fixarray;

using ib::FieldVector;
using ib::Frame;
using ib::Position;
using ib::Result;
using ib::Rotations;
using ib::Status;

namespace {

// ---- goldens harvested from the reference ------------------------------------------------------

/// One `get_mlt1_(iyr, idoy, secs, xGEO, MLT)` observation.
struct MltGolden {
    int year;               ///< `iyr`.
    int doy;                ///< `idoy`.
    double ut;              ///< `secs`.
    fx::vec3d geo;          ///< `xGEO`, Earth radii.
    double mlt;             ///< `MLT`, hours, as the reference returned it.
};

/// One `coord_trans_vec1_(1, sysaxesIN, sysaxesOUT, iyr, idoy, secs, xIN, xOUT)` observation.
struct CoordGolden {
    int year;               ///< `iyr`.
    int doy;                ///< `idoy`.
    double ut;              ///< `secs`.
    int sysaxes_in;         ///< `sysaxesIN`.
    int sysaxes_out;        ///< `sysaxesOUT`.
    fx::vec3d in;           ///< `xIN`, in the input frame's units.
    fx::vec3d out;          ///< `xOUT`, as the reference returned it.
};

// MLT goldens: get_mlt1_(iyr,idoy,secs,xGEO,MLT) on libirbem-O2.so
constexpr MltGolden kMltGoldens[] = {
    {1965, 1, 0.0, fx::vec3d{6.6, 0, 0}, 23.630606468980},
    {1965, 1, 0.0, fx::vec3d{0, 6.6, 0}, 5.580930424733},
    {1965, 1, 0.0, fx::vec3d{-4.2, 1.1, 3.3}, 10.223342693954},
    {1965, 1, 0.0, fx::vec3d{1.5, -1.5, 4}, 21.565814682181},
    {1965, 1, 0.0, fx::vec3d{2, 3, -5}, 2.319914683980},
    {1980, 80, 3600.0, fx::vec3d{6.6, 0, 0}, 0.861424290228},
    {1980, 80, 3600.0, fx::vec3d{0, 6.6, 0}, 6.815930323144},
    {1980, 80, 3600.0, fx::vec3d{-4.2, 1.1, 3.3}, 11.459066258082},
    {1980, 80, 3600.0, fx::vec3d{1.5, -1.5, 4}, 22.788283143732},
    {1980, 80, 3600.0, fx::vec3d{2, 3, -5}, 3.588924890016},
    {1996, 241, 60360.0, fx::vec3d{6.6, 0, 0}, 16.767965034116},
    {1996, 241, 60360.0, fx::vec3d{0, 6.6, 0}, 22.728246130364},
    {1996, 241, 60360.0, fx::vec3d{-4.2, 1.1, 3.3}, 3.381689305107},
    {1996, 241, 60360.0, fx::vec3d{1.5, -1.5, 4}, 14.644901850419},
    {1996, 241, 60360.0, fx::vec3d{2, 3, -5}, 19.562919964898},
    {2003, 300, 70000.0, fx::vec3d{6.6, 0, 0}, 19.655573338726},
    {2003, 300, 70000.0, fx::vec3d{0, 6.6, 0}, 1.618535893965},
    {2003, 300, 70000.0, fx::vec3d{-4.2, 1.1, 3.3}, 6.278388711018},
    {2003, 300, 70000.0, fx::vec3d{1.5, -1.5, 4}, 17.504177503200},
    {2003, 300, 70000.0, fx::vec3d{2, 3, -5}, 22.484768398660},
    {2015, 180, 43200.0, fx::vec3d{6.6, 0, 0}, 11.674330368781},
    {2015, 180, 43200.0, fx::vec3d{0, 6.6, 0}, 17.643231682108},
    {2015, 180, 43200.0, fx::vec3d{-4.2, 1.1, 3.3}, 22.316306992612},
    {2015, 180, 43200.0, fx::vec3d{1.5, -1.5, 4}, 9.468671598061},
    {2015, 180, 43200.0, fx::vec3d{2, 3, -5}, 14.579943533910},
    {2020, 90, 21600.0, fx::vec3d{6.6, 0, 0}, 5.942969050066},
    {2020, 90, 21600.0, fx::vec3d{0, 6.6, 0}, 11.913671666583},
    {2020, 90, 21600.0, fx::vec3d{-4.2, 1.1, 3.3}, 16.594262388902},
    {2020, 90, 21600.0, fx::vec3d{1.5, -1.5, 4}, 3.708466867771},
    {2020, 90, 21600.0, fx::vec3d{2, 3, -5}, 8.876866843900},
    {2024, 355, 86399.0, fx::vec3d{6.6, 0, 0}, 23.773069810563},
    {2024, 355, 86399.0, fx::vec3d{0, 6.6, 0}, 5.744874711923},
    {2024, 355, 86399.0, fx::vec3d{-4.2, 1.1, 3.3}, 10.429919004722},
    {2024, 355, 86399.0, fx::vec3d{1.5, -1.5, 4}, 21.521982999621},
    {2024, 355, 86399.0, fx::vec3d{2, 3, -5}, 2.724197247064},
    {2029, 172, 7200.0, fx::vec3d{6.6, 0, 0}, 2.131654471894},
    {2029, 172, 7200.0, fx::vec3d{0, 6.6, 0}, 8.104952597371},
    {2029, 172, 7200.0, fx::vec3d{-4.2, 1.1, 3.3}, 12.795158808323},
    {2029, 172, 7200.0, fx::vec3d{1.5, -1.5, 4}, 23.861989969113},
    {2029, 172, 7200.0, fx::vec3d{2, 3, -5}, 5.105194445541},
};

// coord_trans goldens: coord_trans_vec1_(1,sysin,sysout,iyr,idoy,secs,xIN,xOUT)
constexpr CoordGolden kCoordGoldens[] = {
    {1965, 1, 0.0, 1, 0, fx::vec3d{-4.2, 1.1, 3.3}, fx::vec3d{28374.6436269, 37.2716465677, 165.323606863}},
    {1965, 1, 0.0, 1, 1, fx::vec3d{-4.2, 1.1, 3.3}, fx::vec3d{-4.2, 1.1, 3.3}},
    {1965, 1, 0.0, 1, 2, fx::vec3d{-4.2, 1.1, 3.3}, fx::vec3d{2.55821034102, -2.11437724681, 4.32723568913}},
    {1965, 1, 0.0, 1, 3, fx::vec3d{-4.2, 1.1, 3.3}, fx::vec3d{2.55821034102, -0.794595569011, 4.75013668371}},
    {1965, 1, 0.0, 1, 4, fx::vec3d{-4.2, 1.1, 3.3}, fx::vec3d{4.21315231194, -2.11437724681, 2.74203509361}},
    {1965, 1, 0.0, 1, 5, fx::vec3d{-4.2, 1.1, 3.3}, fx::vec3d{-0.320937383753, -4.32978050202, 3.3}},
    {1965, 1, 0.0, 1, 6, fx::vec3d{-4.2, 1.1, 3.3}, fx::vec3d{-3.08350807474, -3.56556047466, 2.74203509361}},
    {1965, 1, 0.0, 1, 7, fx::vec3d{-4.2, 1.1, 3.3}, fx::vec3d{5.45343928177, 37.2376730778, 165.323606863}},
    {1965, 1, 0.0, 1, 8, fx::vec3d{-4.2, 1.1, 3.3}, fx::vec3d{5.45343928177, 37.2716465677, 165.323606863}},
    {1980, 80, 3600.0, 1, 0, fx::vec3d{-4.2, 1.1, 3.3}, fx::vec3d{28374.6436269, 37.2716465677, 165.323606863}},
    {1980, 80, 3600.0, 1, 1, fx::vec3d{-4.2, 1.1, 3.3}, fx::vec3d{-4.2, 1.1, 3.3}},
    {1980, 80, 3600.0, 1, 2, fx::vec3d{-4.2, 1.1, 3.3}, fx::vec3d{4.33042648942, -0.663200867145, 3.24770242319}},
    {1980, 80, 3600.0, 1, 3, fx::vec3d{-4.2, 1.1, 3.3}, fx::vec3d{4.33042648942, 1.20841526112, 3.08669634084}},
    {1980, 80, 3600.0, 1, 4, fx::vec3d{-4.2, 1.1, 3.3}, fx::vec3d{4.65174014757, -0.663200867145, 2.7679375371}},
    {1980, 80, 3600.0, 1, 5, fx::vec3d{-4.2, 1.1, 3.3}, fx::vec3d{4.33912922562, -0.148180846763, 3.3}},
    {1980, 80, 3600.0, 1, 6, fx::vec3d{-4.2, 1.1, 3.3}, fx::vec3d{-3.0155488467, -3.60346873774, 2.7679375371}},
    {1980, 80, 3600.0, 1, 7, fx::vec3d{-4.2, 1.1, 3.3}, fx::vec3d{5.45343928177, 37.2376730778, 165.323606863}},
    {1980, 80, 3600.0, 1, 8, fx::vec3d{-4.2, 1.1, 3.3}, fx::vec3d{5.45343928177, 37.2716465677, 165.323606863}},
    {1996, 241, 60360.0, 1, 0, fx::vec3d{-4.2, 1.1, 3.3}, fx::vec3d{28374.6436269, 37.2716465677, 165.323606863}},
    {1996, 241, 60360.0, 1, 1, fx::vec3d{-4.2, 1.1, 3.3}, fx::vec3d{-4.2, 1.1, 3.3}},
    {1996, 241, 60360.0, 1, 2, fx::vec3d{-4.2, 1.1, 3.3}, fx::vec3d{-1.81885130691, -3.62006768682, 3.650601302}},
    {1996, 241, 60360.0, 1, 3, fx::vec3d{-4.2, 1.1, 3.3}, fx::vec3d{-1.81885130691, -4.70976279647, 2.06186438383}},
    {1996, 241, 60360.0, 1, 4, fx::vec3d{-4.2, 1.1, 3.3}, fx::vec3d{-2.96031847102, -3.62006768682, 2.80564154748}},
    {1996, 241, 60360.0, 1, 5, fx::vec3d{-4.2, 1.1, 3.3}, fx::vec3d{3.59923377454, 2.42806841671, 3.3}},
    {1996, 241, 60360.0, 1, 6, fx::vec3d{-4.2, 1.1, 3.3}, fx::vec3d{-2.94526041448, -3.63232936254, 2.80564154748}},
    {1996, 241, 60360.0, 1, 7, fx::vec3d{-4.2, 1.1, 3.3}, fx::vec3d{5.45343928177, 37.2376730778, 165.323606863}},
    {1996, 241, 60360.0, 1, 8, fx::vec3d{-4.2, 1.1, 3.3}, fx::vec3d{5.45343928177, 37.2716465677, 165.323606863}},
    {2003, 300, 70000.0, 1, 0, fx::vec3d{-4.2, 1.1, 3.3}, fx::vec3d{28374.6436269, 37.2716465677, 165.323606863}},
    {2003, 300, 70000.0, 1, 1, fx::vec3d{-4.2, 1.1, 3.3}, fx::vec3d{-4.2, 1.1, 3.3}},
    {2003, 300, 70000.0, 1, 2, fx::vec3d{-4.2, 1.1, 3.3}, fx::vec3d{0.0756604970654, -4.65292002249, 2.84334499375}},
    {2003, 300, 70000.0, 1, 3, fx::vec3d{-4.2, 1.1, 3.3}, fx::vec3d{0.0756604970654, -5.16150541915, 1.75916451706}},
    {2003, 300, 70000.0, 1, 4, fx::vec3d{-4.2, 1.1, 3.3}, fx::vec3d{0.339715801914, -4.65292002249, 2.82399157899}},
    {2003, 300, 70000.0, 1, 5, fx::vec3d{-4.2, 1.1, 3.3}, fx::vec3d{-2.94622783069, 3.1890032251, 3.3}},
    {2003, 300, 70000.0, 1, 6, fx::vec3d{-4.2, 1.1, 3.3}, fx::vec3d{-2.91366054052, -3.64357706332, 2.82399157899}},
    {2003, 300, 70000.0, 1, 7, fx::vec3d{-4.2, 1.1, 3.3}, fx::vec3d{5.45343928177, 37.2376730778, 165.323606863}},
    {2003, 300, 70000.0, 1, 8, fx::vec3d{-4.2, 1.1, 3.3}, fx::vec3d{5.45343928177, 37.2716465677, 165.323606863}},
    {2015, 180, 43200.0, 1, 0, fx::vec3d{-4.2, 1.1, 3.3}, fx::vec3d{28374.6436269, 37.2716465677, 165.323606863}},
    {2015, 180, 43200.0, 1, 1, fx::vec3d{-4.2, 1.1, 3.3}, fx::vec3d{-4.2, 1.1, 3.3}},
    {2015, 180, 43200.0, 1, 2, fx::vec3d{-4.2, 1.1, 3.3}, fx::vec3d{-2.54255543943, 1.97936578243, 4.39971850654}},
    {2015, 180, 43200.0, 1, 3, fx::vec3d{-4.2, 1.1, 3.3}, fx::vec3d{-2.54255543943, 0.898073256984, 4.74038419873}},
    {2015, 180, 43200.0, 1, 4, fx::vec3d{-4.2, 1.1, 3.3}, fx::vec3d{-4.19583188716, 1.97936578243, 2.86654947176}},
    {2015, 180, 43200.0, 1, 5, fx::vec3d{-4.2, 1.1, 3.3}, fx::vec3d{-0.560979504322, -4.30526445131, 3.3}},
    {2015, 180, 43200.0, 1, 6, fx::vec3d{-4.2, 1.1, 3.3}, fx::vec3d{-2.82544256575, -3.67964243285, 2.86654947176}},
    {2015, 180, 43200.0, 1, 7, fx::vec3d{-4.2, 1.1, 3.3}, fx::vec3d{5.45343928177, 37.2376730778, 165.323606863}},
    {2015, 180, 43200.0, 1, 8, fx::vec3d{-4.2, 1.1, 3.3}, fx::vec3d{5.45343928177, 37.2716465677, 165.323606863}},
    {2020, 90, 21600.0, 1, 0, fx::vec3d{-4.2, 1.1, 3.3}, fx::vec3d{28374.6436269, 37.2716465677, 165.323606863}},
    {2020, 90, 21600.0, 1, 1, fx::vec3d{-4.2, 1.1, 3.3}, fx::vec3d{-4.2, 1.1, 3.3}},
    {2020, 90, 21600.0, 1, 2, fx::vec3d{-4.2, 1.1, 3.3}, fx::vec3d{1.4054550335, 4.32063590653, 3.01609043499}},
    {2020, 90, 21600.0, 1, 3, fx::vec3d{-4.2, 1.1, 3.3}, fx::vec3d{1.4054550335, 5.10339435966, 1.31203742616}},
    {2020, 90, 21600.0, 1, 4, fx::vec3d{-4.2, 1.1, 3.3}, fx::vec3d{1.66598739298, 4.32063590653, 2.880380421}},
    {2020, 90, 21600.0, 1, 5, fx::vec3d{-4.2, 1.1, 3.3}, fx::vec3d{0.498001552309, 4.31300295083, 3.3}},
    {2020, 90, 21600.0, 1, 6, fx::vec3d{-4.2, 1.1, 3.3}, fx::vec3d{-2.80784405039, -3.68231183077, 2.880380421}},
    {2020, 90, 21600.0, 1, 7, fx::vec3d{-4.2, 1.1, 3.3}, fx::vec3d{5.45343928177, 37.2376730778, 165.323606863}},
    {2020, 90, 21600.0, 1, 8, fx::vec3d{-4.2, 1.1, 3.3}, fx::vec3d{5.45343928177, 37.2716465677, 165.323606863}},
    {2024, 355, 86399.0, 1, 0, fx::vec3d{-4.2, 1.1, 3.3}, fx::vec3d{28374.6436269, 37.2716465677, 165.323606863}},
    {2024, 355, 86399.0, 1, 1, fx::vec3d{-4.2, 1.1, 3.3}, fx::vec3d{-4.2, 1.1, 3.3}},
    {2024, 355, 86399.0, 1, 2, fx::vec3d{-4.2, 1.1, 3.3}, fx::vec3d{2.54940537475, -1.84813732651, 4.45251846235}},
    {2024, 355, 86399.0, 1, 3, fx::vec3d{-4.2, 1.1, 3.3}, fx::vec3d{2.54940537475, -1.07873776843, 4.6989083718}},
    {2024, 355, 86399.0, 1, 4, fx::vec3d{-4.2, 1.1, 3.3}, fx::vec3d{4.24005633004, -1.84813732651, 2.88899822437}},
    {2024, 355, 86399.0, 1, 5, fx::vec3d{-4.2, 1.1, 3.3}, fx::vec3d{-1.09616920355, -4.20100143742, 3.3}},
    {2024, 355, 86399.0, 1, 6, fx::vec3d{-4.2, 1.1, 3.3}, fx::vec3d{-2.79536757741, -3.68505215251, 2.88899822437}},
    {2024, 355, 86399.0, 1, 7, fx::vec3d{-4.2, 1.1, 3.3}, fx::vec3d{5.45343928177, 37.2376730778, 165.323606863}},
    {2024, 355, 86399.0, 1, 8, fx::vec3d{-4.2, 1.1, 3.3}, fx::vec3d{5.45343928177, 37.2716465677, 165.323606863}},
    {2029, 172, 7200.0, 1, 0, fx::vec3d{-4.2, 1.1, 3.3}, fx::vec3d{28374.6436269, 37.2716465677, 165.323606863}},
    {2029, 172, 7200.0, 1, 1, fx::vec3d{-4.2, 1.1, 3.3}, fx::vec3d{-4.2, 1.1, 3.3}},
    {2029, 172, 7200.0, 1, 2, fx::vec3d{-4.2, 1.1, 3.3}, fx::vec3d{5.16250780898, 0.954395739928, 1.47568353443}},
    {2029, 172, 7200.0, 1, 3, fx::vec3d{-4.2, 1.1, 3.3}, fx::vec3d{5.16250780898, 1.11438135065, 1.35856131882}},
    {2029, 172, 7200.0, 1, 4, fx::vec3d{-4.2, 1.1, 3.3}, fx::vec3d{4.51822963125, 0.954395739928, 2.90081536311}},
    {2029, 172, 7200.0, 1, 5, fx::vec3d{-4.2, 1.1, 3.3}, fx::vec3d{-1.11497290712, 4.19604997782, 3.3}},
    {2029, 172, 7200.0, 1, 6, fx::vec3d{-4.2, 1.1, 3.3}, fx::vec3d{-2.77329117721, -3.69244177687, 2.90081536311}},
    {2029, 172, 7200.0, 1, 7, fx::vec3d{-4.2, 1.1, 3.3}, fx::vec3d{5.45343928177, 37.2376730778, 165.323606863}},
    {2029, 172, 7200.0, 1, 8, fx::vec3d{-4.2, 1.1, 3.3}, fx::vec3d{5.45343928177, 37.2716465677, 165.323606863}},
    {2015, 180, 43200.0, 0, 1, fx::vec3d{500, 33, -75}, fx::vec3d{0.23455110641, -0.875356646093, 0.584862845062}},
    {2015, 180, 43200.0, 0, 2, fx::vec3d{500, 33, -75}, fx::vec3d{0.434171788285, -0.784345887424, 0.599686930031}},
    {2015, 180, 43200.0, 0, 4, fx::vec3d{500, 33, -75}, fx::vec3d{0.131972528818, -0.784345887424, 0.728500382587}},
    {2015, 180, 43200.0, 7, 1, fx::vec3d{1.3, 33, -75}, fx::vec3d{0.282183090217, -1.05312162973, 0.70803074552}},
    {2015, 180, 43200.0, 7, 2, fx::vec3d{1.3, 33, -75}, fx::vec3d{0.524075758281, -0.942907938161, 0.725444153422}},
    {2015, 180, 43200.0, 7, 4, fx::vec3d{1.3, 33, -75}, fx::vec3d{0.158616476301, -0.942907938161, 0.880775472864}},
    {2015, 180, 43200.0, 8, 1, fx::vec3d{1.3, 33, -75}, fx::vec3d{0.282614242706, -1.05473071272, 0.70545893818}},
    {2015, 180, 43200.0, 8, 2, fx::vec3d{1.3, 33, -75}, fx::vec3d{0.523435551683, -0.944947689298, 0.723248980452}},
    {2015, 180, 43200.0, 8, 4, fx::vec3d{1.3, 33, -75}, fx::vec3d{0.158989049593, -0.944947689298, 0.878519405933}},
};

// ---- helpers -----------------------------------------------------------------------------------

/// The IGRF model for a (year, doy, ut) epoch — the argument @ref api::rotations_at needs.
/// Aborts on an epoch outside IGRF-14's window, which is a fixture bug rather than a result.
[[nodiscard]] ib::Igrf<> model_at(int year, int doy, double ut) {
    const ib::DateTime dt = api::doy_and_ut2date_and_time(year, doy, ut);
    const std::optional<ib::Igrf<>> m = ib::Igrf<>::at(
        api::date_and_time2decy(dt.year, dt.month, dt.day, dt.hour, dt.minute, dt.second));
    EXPECT_TRUE(m.has_value());
    return m.value_or(ib::Igrf<>{});
}

/// The rotations for a (year, doy, ut) epoch, built through the adapter's own entry point.
[[nodiscard]] Rotations rotations_for(int year, int doy, double ut) {
    const ib::Igrf<> m = model_at(year, doy, ut);
    const Result<Rotations> r = api::rotations_at(year, doy, ut, m);
    EXPECT_EQ(r.status, Status::Ok);
    return r.value;
}

/// The reference epoch used wherever one epoch will do: 2015-06-29 12:00 UT.
constexpr int kYear = 2015;
constexpr int kDoy = 180;
constexpr double kUt = 43200.0;

/// The nine `sysaxes` codes @ref api::coord_trans accepts.
constexpr std::array<int, 9> kSysaxes{0, 1, 2, 3, 4, 5, 6, 7, 8};

/// Largest componentwise gap between two triples, scaled by the larger component so that a
/// kilometre of GDZ altitude and an Earth radius of GEO are judged on the same footing.
[[nodiscard]] double relative_gap(const fx::vec3d& a, const fx::vec3d& b) {
    double worst = 0.0;
    for (std::size_t i = 0; i < 3; ++i) {
        const double scale = std::max({std::abs(a[i]), std::abs(b[i]), 1.0});
        worst = std::max(worst, std::abs(a[i] - b[i]) / scale);
    }
    return worst;
}

}  // namespace

// ---- library information -----------------------------------------------------------------------

TEST(IrbemApi, LibraryInfoMatchesTheReference) {
    // Measured on libirbem-O2.so: get_irbem_ntime_max1_ -> 100000, get_igrf_version_ -> 14.
    EXPECT_EQ(api::get_irbem_ntime_max(), 100000);
    EXPECT_EQ(api::get_igrf_version(), 14);

    // These two deliberately do NOT match the reference, and that is the assertion. The reference
    // returns its own Fortran revision (measured: version 1, release "e7cecb0" space-padded to 80
    // characters); reporting a plausible-looking IRBEM revision from a C++ implementation would let
    // a caller feature-detect against a number that means nothing about what is running.
    EXPECT_EQ(api::irbem_fortran_version(), 1);
    EXPECT_NE(api::irbem_fortran_release(), "e7cecb0");
    EXPECT_NE(api::irbem_fortran_release().find("not IRBEM Fortran"), std::string_view::npos);

    // One definition, two spellings: the adapter must forward, not re-type the literal.
    EXPECT_EQ(api::get_irbem_ntime_max(), ib::max_batch_times());
    EXPECT_EQ(api::get_igrf_version(), ib::igrf_generation());
    EXPECT_EQ(api::irbem_fortran_version(), ib::implementation_version());
    EXPECT_EQ(api::irbem_fortran_release(), ib::implementation_release());
}

// ---- the epoch object ---------------------------------------------------------------------------

TEST(IrbemApi, RotationsAreBuiltOncePerEpoch) {
    const ib::Igrf<> m = model_at(kYear, kDoy, kUt);
    const Result<Rotations> a = api::rotations_at(kYear, kDoy, kUt, m);
    const Result<Rotations> b = api::rotations_at(kYear, kDoy, kUt, m);
    ASSERT_EQ(a.status, Status::Ok);
    ASSERT_EQ(b.status, Status::Ok);

    // Bitwise identical, not merely close: the whole point of hoisting is that the SAME matrices
    // are reused, so a second build must reproduce the first exactly or the hoist changes answers.
    EXPECT_EQ(a.value.jd_ut1, b.value.jd_ut1);
    EXPECT_EQ(a.value.gmst_deg, b.value.gmst_deg);
    EXPECT_EQ(a.value.dipole_tilt_deg, b.value.dipole_tilt_deg);
    for (std::size_t c = 0; c < 3; ++c) {
        for (std::size_t r = 0; r < 3; ++r) {
            EXPECT_EQ(a.value.geo_to_gsm(r, c), b.value.geo_to_gsm(r, c));
        }
    }

    // An epoch IGRF-14 does not cover is reported, not clamped to the nearest year it does.
    const Result<Rotations> early = api::rotations_at(1850, 1, 0.0, m);
    EXPECT_EQ(early.status, Status::DomainError);
    const Result<Rotations> late = api::rotations_at(2099, 1, 0.0, m);
    EXPECT_EQ(late.status, Status::DomainError);
}

// ---- the named transforms ------------------------------------------------------------------------

TEST(IrbemApi, TransformsRoundTrip) {
    const Rotations r = rotations_for(kYear, kDoy, kUt);
    const Position<Frame::GEO> geo{fx::vec3d{-4.2, 1.1, 3.3}};

    // Every named pair, forward then back. A frame rotation is orthogonal, so the inverse is the
    // transpose and the round trip is exact to a rounding or two rather than to a tolerance that
    // hides a wrong matrix.
    constexpr double kExact = 1e-15;
    EXPECT_LT(relative_gap(api::gsm2geo(api::geo2gsm(geo, r), r).v, geo.v), kExact);
    EXPECT_LT(relative_gap(api::sm2geo(api::geo2sm(geo, r), r).v, geo.v), kExact);
    EXPECT_LT(relative_gap(api::gse2geo(api::geo2gse(geo, r), r).v, geo.v), kExact);
    EXPECT_LT(relative_gap(api::mag2geo(api::geo2mag(geo, r), r).v, geo.v), kExact);
    EXPECT_LT(relative_gap(api::gei2geo(api::geo2gei(geo, r), r).v, geo.v), kExact);

    const Position<Frame::GSM> gsm = api::geo2gsm(geo, r);
    EXPECT_LT(relative_gap(api::sm2gsm(api::gsm2sm(gsm, r), r).v, gsm.v), kExact);

    // And each name is exactly its typed counterpart — the adapter adds a spelling, not a formula.
    EXPECT_EQ(api::geo2gsm(geo, r), ib::transform<Frame::GSM>(geo, r));
    EXPECT_EQ(api::geo2sm(geo, r), ib::transform<Frame::SM>(geo, r));
    EXPECT_EQ(api::geo2gse(geo, r), ib::transform<Frame::GSE>(geo, r));
    EXPECT_EQ(api::geo2mag(geo, r), ib::transform<Frame::MAG>(geo, r));
    EXPECT_EQ(api::geo2gei(geo, r), ib::transform<Frame::GEI>(geo, r));
    EXPECT_EQ(api::gsm2sm(gsm, r), ib::transform<Frame::SM>(gsm, r));
    EXPECT_EQ(api::gsm2geo(gsm, r), ib::transform<Frame::GEO>(gsm, r));
    EXPECT_EQ(api::sm2gsm(api::gsm2sm(gsm, r), r), ib::transform<Frame::GSM>(api::gsm2sm(gsm, r), r));
    EXPECT_EQ(api::sm2geo(api::geo2sm(geo, r), r), ib::transform<Frame::GEO>(api::geo2sm(geo, r), r));
    EXPECT_EQ(api::gse2geo(api::geo2gse(geo, r), r), ib::transform<Frame::GEO>(api::geo2gse(geo, r), r));
    EXPECT_EQ(api::mag2geo(api::geo2mag(geo, r), r), ib::transform<Frame::GEO>(api::geo2mag(geo, r), r));
    EXPECT_EQ(api::gei2geo(api::geo2gei(geo, r), r), ib::transform<Frame::GEO>(api::geo2gei(geo, r), r));
}

TEST(IrbemApi, GeodeticRoundTrips) {
    // The equator on the ellipsoid: altitude 0, latitude 0, longitude 0 puts the point on +X at
    // exactly the semi-major axis, so the assertion can be `==` rather than approximate.
    const Position<Frame::GDZ> on_equator{fx::vec3d{0.0, 0.0, 0.0}};
    const Position<Frame::GEO> geo = api::gdz2geo(on_equator);
    EXPECT_EQ(geo.v[1], 0.0);
    EXPECT_EQ(geo.v[2], 0.0);
    EXPECT_EQ(geo.v[0], ib::detail::wgs84_semi_major_km / ib::detail::earth_radius_km);

    for (double alt : {0.0, 500.0, 35786.0}) {
        for (double lat : {-89.0, -45.0, 0.0, 12.5, 45.0, 89.0}) {
            const Position<Frame::GDZ> gdz{fx::vec3d{alt, lat, -75.0}};
            const Position<Frame::GDZ> back = api::geo2gdz(api::gdz2geo(gdz));
            EXPECT_LT(relative_gap(back.v, gdz.v), 1e-12) << alt << " " << lat;

            // RLL is GDZ re-expressed with geocentric radius, so rll2gdz must undo that exactly.
            const Position<Frame::GEO> cart = api::gdz2geo(gdz);
            const Position<Frame::RLL> rll{fx::vec3d{cart.radius(), gdz.latitude(), gdz.longitude()}};
            EXPECT_NEAR(api::rll2gdz(rll).radius(), alt, 1e-6) << alt << " " << lat;
        }
    }
}

TEST(IrbemApi, SphericalRoundTrips) {
    // (1, 0, 0) in SPH is the +X unit vector exactly: cos(0)=1 and sin(0)=0 are both exact.
    const Position<Frame::GEO> unit = api::sph2car(Position<Frame::SPH>{fx::vec3d{1.0, 0.0, 0.0}});
    EXPECT_EQ(unit.v[0], 1.0);
    EXPECT_EQ(unit.v[1], 0.0);
    EXPECT_EQ(unit.v[2], 0.0);

    for (double radius : {1.0, 2.5, 6.6}) {
        for (double lat : {-80.0, -30.0, 0.0, 30.0, 80.0}) {
            for (double lon : {-179.0, -75.0, 0.0, 45.0, 179.0}) {
                const Position<Frame::SPH> sph{fx::vec3d{radius, lat, lon}};
                EXPECT_LT(relative_gap(api::car2sph(api::sph2car(sph)).v, sph.v), 1e-13);
            }
        }
    }
}

// ---- the heliospheric pairs ----------------------------------------------------------------------

TEST(IrbemApi, HeliosphericPairsRoundTrip) {
    // The Fränz & Harper reference epoch, 1996-08-28 16:46:00 TT (MJD 50323.69861111).
    const ib::HelioGeometry g = ib::helio_geometry(50323.0 + (60360.0 / 86400.0));
    const fx::vec3d seed{-5.7864918, -3.0028771, 3.3908764};

    // Both vector kinds, because clang counts template instantiations separately and because the
    // rotational pairs are written once over both — a per-kind body is exactly what must not exist.
    const Position<Frame::HAE> hae_p{seed};
    const FieldVector<Frame::HAE> hae_f{seed};

    EXPECT_LT(relative_gap(api::hee2hae(api::hae2hee(hae_p, g), g).v, seed), 1e-15);
    EXPECT_LT(relative_gap(api::hee2hae(api::hae2hee(hae_f, g), g).v, seed), 1e-15);
    EXPECT_LT(relative_gap(api::heeq2hae(api::hae2heeq(hae_p, g), g).v, seed), 1e-15);
    EXPECT_LT(relative_gap(api::heeq2hae(api::hae2heeq(hae_f, g), g).v, seed), 1e-15);

    // A rotation preserves length; an origin shift does not. These four are rotations.
    EXPECT_NEAR(fx::norm(api::hae2hee(hae_p, g).v), fx::norm(seed), 1e-13);
    EXPECT_NEAR(fx::norm(api::hae2heeq(hae_f, g).v), fx::norm(seed), 1e-13);

    // And each name is exactly its canonical counterpart in coords_helio.hpp.
    EXPECT_EQ(api::hae2hee(hae_p, g), ib::HAE2HEE<Position>(hae_p, g));
    EXPECT_EQ(api::hae2heeq(hae_f, g), ib::HAE2HEEQ<FieldVector>(hae_f, g));
    EXPECT_EQ(api::hee2hae(api::hae2hee(hae_p, g), g), ib::HEE2HAE<Position>(api::hae2hee(hae_p, g), g));
    EXPECT_EQ(api::heeq2hae(api::hae2heeq(hae_f, g), g), ib::HEEQ2HAE<FieldVector>(api::hae2heeq(hae_f, g), g));
}

TEST(IrbemApi, HeliosphericOriginShiftAppliesToPositionsOnly) {
    const ib::HelioGeometry g = ib::helio_geometry(50323.0 + (60360.0 / 86400.0));

    // The Earth sits at the GSE origin, and in HEE that is (r0, 0, 0). Exact: the transform is two
    // sign flips and one add, and 0.0 negated is still 0.0 in magnitude.
    const Position<Frame::HEE> earth = api::gse2hee(Position<Frame::GSE>{fx::vec3d{0.0, 0.0, 0.0}}, g);
    EXPECT_EQ(earth.v[0], g.sun_earth_distance_re);
    EXPECT_EQ(earth.v[1], 0.0);
    EXPECT_EQ(earth.v[2], 0.0);

    // A FIELD measured at the Earth must not move at all: it has a direction, not a place. This is
    // the ~23 000 Re defect the two overloads exist to make impossible.
    const FieldVector<Frame::HEE> b = api::gse2hee(FieldVector<Frame::GSE>{fx::vec3d{5.0, -3.0, 2.0}}, g);
    EXPECT_EQ(b.v[0], -5.0);
    EXPECT_EQ(b.v[1], 3.0);
    EXPECT_EQ(b.v[2], 2.0);
    EXPECT_EQ(fx::norm(b.v), fx::norm(fx::vec3d{5.0, -3.0, 2.0}));

    // Both directions, both kinds — the map is an involution for a position and a half turn for a
    // field, so each is its own inverse and the assertion can be `==`.
    EXPECT_EQ(api::hee2gse(earth, g), (Position<Frame::GSE>{fx::vec3d{0.0, 0.0, 0.0}}));
    EXPECT_EQ(api::hee2gse(b, g), (FieldVector<Frame::GSE>{fx::vec3d{5.0, -3.0, 2.0}}));

    // The AU bridge the file brief promises: IRBEM's gse2hee1_ returns 1.01655079 for GSE (1,0,0)
    // Re at 2015-180 12:00 UT. That is this module's answer divided by au_in_earth_radii.
    const ib::HelioGeometry g2015 = ib::helio_geometry(57202.0 + 0.5);
    const Position<Frame::HEE> one = api::gse2hee(Position<Frame::GSE>{fx::vec3d{1.0, 0.0, 0.0}}, g2015);
    EXPECT_NEAR(one.v[0] / ib::au_in_earth_radii(), 1.01655079, 1e-4);
}

// ---- magnetic local time --------------------------------------------------------------------------

TEST(IrbemApi, MltIsTheSolarMagneticClockAngle) {
    // A Rotations built as an aggregate with an identity GEO->SM block makes the clock angle
    // testable exactly: the GEO components ARE the SM components, so the four cardinal azimuths
    // land on 12, 18, 0 and 6 hours with no rounding from an intervening rotation.
    Rotations id{};
    id.geo_to_sm = fx::mat3d::identity();

    EXPECT_EQ(api::get_mlt(Position<Frame::GEO>{fx::vec3d{7.0, 0.0, 1.0}}, id).value, 12.0);
    EXPECT_EQ(api::get_mlt(Position<Frame::GEO>{fx::vec3d{0.0, 7.0, 1.0}}, id).value, 18.0);
    EXPECT_EQ(api::get_mlt(Position<Frame::GEO>{fx::vec3d{0.0, -7.0, 1.0}}, id).value, 6.0);
    // atan2(+0, -x) is +pi, so the raw hour angle is 24 exactly and must fold to 0, not stay at 24.
    const Result<double> midnight = api::get_mlt(Position<Frame::GEO>{fx::vec3d{-7.0, 0.0, 1.0}}, id);
    EXPECT_EQ(midnight.status, Status::Ok);
    EXPECT_EQ(midnight.value, 0.0);
    // ...and from the other side of the tail, atan2(-0, -x) = -pi, which must fold the same way.
    const Result<double> midnight_below =
        api::get_mlt(Position<Frame::GEO>{fx::vec3d{-7.0, -0.0, 1.0}}, id);
    EXPECT_EQ(midnight_below.value, 0.0);

    // Every value is in range, across a full turn.
    for (int i = 0; i < 360; ++i) {
        const double a = static_cast<double>(i) * (std::numbers::pi / 180.0);
        const Result<double> h = api::get_mlt(
            Position<Frame::GEO>{fx::vec3d{4.0 * std::cos(a), 4.0 * std::sin(a), 0.5}}, id);
        ASSERT_EQ(h.status, Status::Ok);
        EXPECT_GE(h.value, 0.0);
        EXPECT_LT(h.value, 24.0);
    }

    // The dipole axis itself has no magnetic meridian, so MLT is undefined there and says so. A
    // default-constructed Rotations sends everything to the SM origin, which is that degenerate
    // locus for every input.
    const Rotations degenerate{};
    const Result<double> on_axis =
        api::get_mlt(Position<Frame::GEO>{fx::vec3d{1.0, 2.0, 3.0}}, degenerate);
    EXPECT_EQ(on_axis.status, Status::DomainError);
    EXPECT_EQ(on_axis.value, 12.0);  // in range, not a NaN, for a caller that ignores the status

    // A non-finite input is reported rather than propagated as a NaN hour.
    const double nan = std::numeric_limits<double>::quiet_NaN();
    const Result<double> bad = api::get_mlt(Position<Frame::GEO>{fx::vec3d{nan, 1.0, 1.0}}, id);
    EXPECT_EQ(bad.status, Status::DomainError);
    EXPECT_EQ(bad.value, 0.0);
}

TEST(IrbemApi, MltMatchesTheReference) {
    // The residual is not arithmetic: the two libraries build the SM frame from different solar
    // ephemeris and sidereal-time series, so they disagree by ~1e-3 h (0.015 deg of SM azimuth)
    // even where both are correct. The cap is set from the observed maximum over the sweep and the
    // maximum is PRINTED, so a regression that doubles it is visible rather than merely under cap.
    double worst = 0.0;
    for (const MltGolden& g : kMltGoldens) {
        const Rotations r = rotations_for(g.year, g.doy, g.ut);
        const Result<double> mine = api::get_mlt(Position<Frame::GEO>{g.geo}, r);
        ASSERT_EQ(mine.status, Status::Ok);
        // Both are hour angles on a 24 h circle, so 23.99 and 0.01 are 0.02 apart, not 23.98.
        double d = std::abs(mine.value - g.mlt);
        if (d > 12.0) d = 24.0 - d;
        worst = std::max(worst, d);
        EXPECT_LT(d, 5e-3) << g.year << "/" << g.doy << " @" << g.ut;
    }
    std::printf("[ MLT vs oracle ] max |dMLT| = %.4g h (%.4g s) over %zu points\n", worst,
                worst * 3600.0, std::size(kMltGoldens));
    EXPECT_GT(worst, 0.0);  // a zero here would mean the goldens were not actually compared
}

TEST(IrbemApi, MltBatchAgreesWithTheScalarLane) {
    std::vector<Position<Frame::GEO>> pts;
    for (int i = 0; i < 64; ++i) {
        const double t = 0.1 * static_cast<double>(i);
        pts.push_back(Position<Frame::GEO>{fx::vec3d{2.0 + t, 1.0 - (0.05 * t), -3.0 + (0.2 * t)}});
    }
    const Rotations r = rotations_for(kYear, kDoy, kUt);
    const std::array<Rotations, 1> shared{r};
    std::vector<double> mlt(pts.size());
    std::vector<Status> st(pts.size());

    const Result<bool> batch = api::get_mlt_vec(pts, shared, mlt, st);
    EXPECT_EQ(batch.status, Status::Ok);
    // The coordinate lanes are host-only by measurement, not by omission — see the file brief.
    EXPECT_FALSE(batch.value);
    for (std::size_t i = 0; i < pts.size(); ++i) {
        // `==`, not NEAR: the batch must compute the identical expression, not an equivalent one.
        EXPECT_EQ(mlt[i], api::get_mlt(pts[i], r).value);
        EXPECT_EQ(st[i], Status::Ok);
    }

    // Per-point epochs: a multi-year survey, where hoisting is not available and must not be faked.
    std::vector<Rotations> per_point;
    for (std::size_t i = 0; i < pts.size(); ++i) {
        per_point.push_back(rotations_for(2000 + static_cast<int>(i % 25), 1 + static_cast<int>(i), 0.0));
    }
    const Result<bool> varied = api::get_mlt_vec(pts, per_point, mlt, st);
    EXPECT_EQ(varied.status, Status::Ok);
    for (std::size_t i = 0; i < pts.size(); ++i) {
        EXPECT_EQ(mlt[i], api::get_mlt(pts[i], per_point[i]).value);
    }

    // A degenerate point is reported per point, and does not poison its neighbours.
    const std::array<Rotations, 1> degenerate{Rotations{}};
    const Result<bool> some_bad = api::get_mlt_vec(pts, degenerate, mlt, st);
    EXPECT_EQ(some_bad.status, Status::DomainError);
    EXPECT_EQ(st[0], Status::DomainError);
}

// ---- the runtime frame-pair dispatcher ------------------------------------------------------------

TEST(IrbemApi, CoordTransHubRejectsHeliosphericFrames) {
    // The heliospheric frames need a HelioGeometry, not a Rotations, so the hub reports them rather
    // than producing a plausible-looking triple. They are unreachable through coord_trans — sysaxes
    // stops at 8 — so the hub is driven directly here, which is the only way this line runs.
    const Rotations r = rotations_for(kYear, kDoy, kUt);
    const fx::vec3d v{1.0, 2.0, 3.0};
    for (Frame f : {Frame::HEE, Frame::HAE, Frame::HEEQ}) {
        EXPECT_EQ(api::detail::to_geo(f, v, r).status, Status::DomainError);
        EXPECT_EQ(api::detail::from_geo(f, v, r).status, Status::DomainError);
        EXPECT_EQ(api::detail::to_geo(f, v, r).value, fx::vec3d{});
        EXPECT_EQ(api::detail::from_geo(f, v, r).value, fx::vec3d{});
        // ...and no sysaxes code names them, which is what makes the above unreachable in use.
        EXPECT_FALSE(ib::sysaxes_of(f).has_value());
    }
}

TEST(IrbemApi, CoordTransCoversEveryFramePair) {
    const Rotations r = rotations_for(kYear, kDoy, kUt);
    const fx::vec3d seed_geo{-4.2, 1.1, 3.3};

    double worst = 0.0;
    for (int a : kSysaxes) {
        // Express the seed point in frame `a` first, so every leg is fed something physical.
        const Result<fx::vec3d> in = api::coord_trans(1, a, seed_geo, r);
        ASSERT_EQ(in.status, Status::Ok) << a;
        for (int b : kSysaxes) {
            const Result<fx::vec3d> there = api::coord_trans(a, b, in.value, r);
            ASSERT_EQ(there.status, Status::Ok) << a << "->" << b;
            const Result<fx::vec3d> back = api::coord_trans(b, a, there.value, r);
            ASSERT_EQ(back.status, Status::Ok) << b << "->" << a;
            const double d = relative_gap(back.value, in.value);
            worst = std::max(worst, d);
            EXPECT_LT(d, 1e-11) << a << " -> " << b << " -> " << a;

            // Same frame in and out is the identity, exactly — not a round trip through GEO.
            if (a == b) EXPECT_EQ(there.value, in.value);
        }
    }
    std::printf("[ coord_trans ] worst 81-pair round-trip gap = %.4g relative\n", worst);
    EXPECT_GT(worst, 0.0);
}

TEST(IrbemApi, CoordTransReportsBadSysaxes) {
    const Rotations r = rotations_for(kYear, kDoy, kUt);
    const fx::vec3d v{1.5, 2.5, 3.5};

    // 9, 10, 11 name frames this module HAS but that need a HelioGeometry; 12, 13, 14 name TOD,
    // J2000 and TEME, which it does not have at all; -1 and 15 are outside IRBEM's own table. All
    // six are reported the same way, and none of them silently becomes GEO.
    for (int bad : {-1, 9, 10, 11, 12, 13, 14, 15, 1000}) {
        EXPECT_EQ(api::coord_trans(bad, 1, v, r).status, Status::DomainError) << bad;
        EXPECT_EQ(api::coord_trans(1, bad, v, r).status, Status::DomainError) << bad;
        EXPECT_EQ(api::coord_trans(bad, 1, v, r).value, fx::vec3d{}) << bad;
    }

    // A non-finite component is a caller error too, and is caught before any arithmetic runs.
    const double inf = std::numeric_limits<double>::infinity();
    const double nan = std::numeric_limits<double>::quiet_NaN();
    EXPECT_EQ(api::coord_trans(1, 2, fx::vec3d{nan, 1.0, 1.0}, r).status, Status::DomainError);
    EXPECT_EQ(api::coord_trans(1, 2, fx::vec3d{1.0, inf, 1.0}, r).status, Status::DomainError);
    EXPECT_EQ(api::coord_trans(1, 2, fx::vec3d{1.0, 1.0, nan}, r).status, Status::DomainError);
}

TEST(IrbemApi, CoordTransMatchesTheReference) {
    // Nine output frames × eight epochs from 1965 to 2029, plus a few non-GEO inputs. The residual
    // is again two implementations of the same geometry, not two geometries: the rotational legs
    // differ by the solar-ephemeris series (~1e-5 relative) and the geodetic legs by the
    // reference's iterative latitude solve against this module's closed form (~2e-7 deg, a fact
    // coords_geodetic.hpp already documents and measures).
    double worst = 0.0;
    int worst_pair[2] = {0, 0};
    for (const CoordGolden& g : kCoordGoldens) {
        const Rotations r = rotations_for(g.year, g.doy, g.ut);
        const Result<fx::vec3d> mine = api::coord_trans(g.sysaxes_in, g.sysaxes_out, g.in, r);
        ASSERT_EQ(mine.status, Status::Ok);
        const double d = relative_gap(mine.value, g.out);
        if (d > worst) {
            worst = d;
            worst_pair[0] = g.sysaxes_in;
            worst_pair[1] = g.sysaxes_out;
        }
        EXPECT_LT(d, 5e-5) << g.year << "/" << g.doy << " " << g.sysaxes_in << "->" << g.sysaxes_out;
    }
    std::printf("[ coord_trans vs oracle ] worst = %.4g relative, at sysaxes %d->%d, over %zu rows\n",
                worst, worst_pair[0], worst_pair[1], std::size(kCoordGoldens));
    EXPECT_GT(worst, 0.0);
}

TEST(IrbemApi, CoordTransVecAgreesWithTheScalarLane) {
    const Rotations r = rotations_for(kYear, kDoy, kUt);
    std::vector<fx::vec3d> in;
    for (int i = 0; i < 100; ++i) {
        const double t = 0.05 * static_cast<double>(i);
        in.push_back(fx::vec3d{2.0 + t, -1.0 + (0.3 * t), 0.5 - (0.1 * t)});
    }
    std::vector<fx::vec3d> out(in.size());
    std::vector<Status> st(in.size());
    const std::array<Rotations, 1> shared{r};

    for (int a : kSysaxes) {
        for (int b : kSysaxes) {
            const Result<bool> batch = api::coord_trans_vec(a, b, in, out, shared, st);
            ASSERT_EQ(batch.status, Status::Ok) << a << "->" << b;
            EXPECT_FALSE(batch.value);
            for (std::size_t i = 0; i < in.size(); ++i) {
                EXPECT_EQ(out[i], api::coord_trans(a, b, in[i], r).value);
                EXPECT_EQ(st[i], Status::Ok);
            }
        }
    }

    // A bad frame code fails the whole call and writes nothing: the fault is in the call, not the
    // data, so leaving the output untouched is more useful than filling it with a sentinel.
    std::fill(out.begin(), out.end(), fx::vec3d{-7.0, -7.0, -7.0});
    EXPECT_EQ(api::coord_trans_vec(1, 12, in, out, shared, st).status, Status::DomainError);
    EXPECT_EQ(api::coord_trans_vec(12, 1, in, out, shared, st).status, Status::DomainError);
    EXPECT_EQ(out[0], (fx::vec3d{-7.0, -7.0, -7.0}));

    // A non-finite point is reported per point and the rest of the batch still converts.
    std::vector<fx::vec3d> mixed = in;
    mixed[3] = fx::vec3d{std::numeric_limits<double>::quiet_NaN(), 1.0, 1.0};
    const Result<bool> partial = api::coord_trans_vec(1, 2, mixed, out, shared, st);
    EXPECT_EQ(partial.status, Status::DomainError);
    EXPECT_EQ(st[3], Status::DomainError);
    EXPECT_EQ(st[4], Status::Ok);
    EXPECT_EQ(out[4], api::coord_trans(1, 2, mixed[4], r).value);
}

TEST(IrbemApi, CoordTransVecTakesPerPointEpochs) {
    // IRBEM's COORD_TRANS_VEC takes a per-point iyr/idoy/secs, and a multi-year survey needs it.
    std::vector<fx::vec3d> in;
    std::vector<Rotations> per_point;
    for (int i = 0; i < 32; ++i) {
        in.push_back(fx::vec3d{3.0 + (0.1 * i), 1.0, -2.0});
        per_point.push_back(rotations_for(1990 + i, 1 + (10 * i), 900.0 * i));
    }
    std::vector<fx::vec3d> out(in.size());
    std::vector<Status> st(in.size());

    const Result<bool> batch = api::coord_trans_vec(1, 4, in, out, per_point, st);
    ASSERT_EQ(batch.status, Status::Ok);
    for (std::size_t i = 0; i < in.size(); ++i) {
        EXPECT_EQ(out[i], api::coord_trans(1, 4, in[i], per_point[i]).value);
    }
    // The epochs really are different, so a shared-epoch implementation would fail this.
    EXPECT_NE(out.front(), out.back());
    EXPECT_NE(per_point.front().gmst_deg, per_point.back().gmst_deg);
}

TEST(IrbemApi, TypedBatchAgreesWithTheScalarTransform) {
    const Rotations r = rotations_for(kYear, kDoy, kUt);
    std::vector<Position<Frame::GEO>> geo;
    std::vector<FieldVector<Frame::GEO>> field;
    for (int i = 0; i < 48; ++i) {
        const double t = 0.25 * static_cast<double>(i);
        geo.push_back(Position<Frame::GEO>{fx::vec3d{1.5 + t, -0.5 + (0.1 * t), 2.0 - (0.05 * t)}});
        field.push_back(FieldVector<Frame::GEO>{fx::vec3d{100.0 * t, -50.0, 25.0 + t}});
    }
    std::vector<Position<Frame::GSM>> gsm(geo.size());
    std::vector<FieldVector<Frame::MAG>> mag(field.size());
    const std::array<Rotations, 1> shared{r};

    // Shared epoch: the matrix is hoisted out of the loop, and hoisting must not change a bit.
    const Result<bool> a = api::transform_vec<Frame::GSM>(
        std::span<const Position<Frame::GEO>>{geo}, std::span<Position<Frame::GSM>>{gsm}, shared);
    ASSERT_EQ(a.status, Status::Ok);
    EXPECT_FALSE(a.value);
    for (std::size_t i = 0; i < geo.size(); ++i) EXPECT_EQ(gsm[i], ib::transform<Frame::GSM>(geo[i], r));

    const Result<bool> b = api::transform_vec<Frame::MAG>(
        std::span<const FieldVector<Frame::GEO>>{field}, std::span<FieldVector<Frame::MAG>>{mag},
        shared);
    ASSERT_EQ(b.status, Status::Ok);
    for (std::size_t i = 0; i < field.size(); ++i) EXPECT_EQ(mag[i], ib::transform<Frame::MAG>(field[i], r));

    // Per-point epochs take the other branch and must agree with the scalar lane just as exactly.
    std::vector<Rotations> per_point;
    for (std::size_t i = 0; i < geo.size(); ++i) {
        per_point.push_back(rotations_for(1980 + static_cast<int>(i), 100, 3600.0));
    }
    std::vector<Position<Frame::GEO>> back(geo.size());
    const Result<bool> c = api::transform_vec<Frame::GEO>(
        std::span<const Position<Frame::GSM>>{gsm}, std::span<Position<Frame::GEO>>{back},
        per_point);
    ASSERT_EQ(c.status, Status::Ok);
    for (std::size_t i = 0; i < gsm.size(); ++i) {
        EXPECT_EQ(back[i], ib::transform<Frame::GEO>(gsm[i], per_point[i]));
    }
}

TEST(IrbemApi, BatchRejectsMismatchedSpans) {
    const Rotations r = rotations_for(kYear, kDoy, kUt);
    const std::array<Rotations, 1> shared{r};
    const std::array<Rotations, 2> two{r, r};

    std::vector<fx::vec3d> in(4, fx::vec3d{1.0, 2.0, 3.0});
    std::vector<fx::vec3d> out(4);
    std::vector<fx::vec3d> short_out(3);
    std::vector<Status> st(4);
    std::vector<Status> short_st(3);
    EXPECT_EQ(api::coord_trans_vec(1, 2, in, short_out, shared, st).status, Status::DomainError);
    EXPECT_EQ(api::coord_trans_vec(1, 2, in, out, shared, short_st).status, Status::DomainError);
    EXPECT_EQ(api::coord_trans_vec(1, 2, in, out, two, st).status, Status::DomainError);

    std::vector<Position<Frame::GEO>> pts(4, Position<Frame::GEO>{fx::vec3d{1.0, 2.0, 3.0}});
    std::vector<double> mlt(4);
    std::vector<double> short_mlt(3);
    EXPECT_EQ(api::get_mlt_vec(pts, shared, short_mlt, st).status, Status::DomainError);
    EXPECT_EQ(api::get_mlt_vec(pts, shared, mlt, short_st).status, Status::DomainError);
    EXPECT_EQ(api::get_mlt_vec(pts, two, mlt, st).status, Status::DomainError);

    std::vector<Position<Frame::GSM>> gsm(4);
    std::vector<Position<Frame::GSM>> short_gsm(3);
    EXPECT_EQ((api::transform_vec<Frame::GSM>(std::span<const Position<Frame::GEO>>{pts},
                                              std::span<Position<Frame::GSM>>{short_gsm}, shared)
                   .status),
              Status::DomainError);
    EXPECT_EQ((api::transform_vec<Frame::GSM>(std::span<const Position<Frame::GEO>>{pts},
                                              std::span<Position<Frame::GSM>>{gsm}, two)
                   .status),
              Status::DomainError);

    // An empty batch is a valid batch: zero points converted, no error.
    std::vector<fx::vec3d> none;
    std::vector<Status> no_status;
    EXPECT_EQ(api::coord_trans_vec(1, 2, none, none, shared, no_status).status, Status::Ok);
}

// ---- date and time -------------------------------------------------------------------------------

TEST(IrbemApi, DateRoundTrips) {
    // A leap year and a common year, at the boundaries where an off-by-one shows up.
    EXPECT_EQ(api::get_doy(2015, 1, 1), 1);
    EXPECT_EQ(api::get_doy(2015, 12, 31), 365);
    EXPECT_EQ(api::get_doy(2016, 12, 31), 366);
    EXPECT_EQ(api::get_doy(2016, 3, 1), 61);   // 31 + 29 + 1
    EXPECT_EQ(api::get_doy(2015, 3, 1), 60);   // 31 + 28 + 1

    // JULDAY/CALDAT are inverses over a long stretch, including the Gregorian-calendar arithmetic
    // either side of a century that is not a leap year.
    for (int year : {1899, 1900, 1901, 1999, 2000, 2001, 2015, 2016, 2029, 2100}) {
        for (int month = 1; month <= 12; ++month) {
            for (int day : {1, 15, 28}) {
                const std::int64_t jdn = api::julday(year, month, day);
                const ib::CalendarDate back = api::caldat(jdn);
                EXPECT_EQ(back.year, year);
                EXPECT_EQ(back.month, month);
                EXPECT_EQ(back.day, day);
            }
        }
    }
    // Consecutive days are consecutive Julian Day Numbers across a month and a year boundary.
    EXPECT_EQ(api::julday(2016, 3, 1) - api::julday(2016, 2, 29), 1);
    EXPECT_EQ(api::julday(2016, 1, 1) - api::julday(2015, 12, 31), 1);

    // The decimal year and its inverse, at an exactly-representable instant: noon on day 1 of a
    // 365-day year is 0.5/365 of the way through it.
    EXPECT_EQ(api::date_and_time2decy(2015, 1, 1, 0, 0, 0), 2015.0);
    EXPECT_NEAR(api::date_and_time2decy(2015, 1, 1, 12, 0, 0), 2015.0 + (0.5 / 365.0), 1e-15);
    const ib::DateTime noon = api::decy2date_and_time(2015.0 + (0.5 / 365.0));
    EXPECT_EQ(noon.year, 2015);
    EXPECT_EQ(noon.month, 1);
    EXPECT_EQ(noon.day, 1);
    EXPECT_EQ(noon.hour, 12);

    // DOY_AND_UT2DATE_AND_TIME closes the loop with GET_DOY.
    for (int doy : {1, 59, 60, 200, 365}) {
        const ib::DateTime dt = api::doy_and_ut2date_and_time(2015, doy, 3661.0);
        EXPECT_EQ(dt.day_of_year, doy);
        EXPECT_EQ(api::get_doy(dt.year, dt.month, dt.day), doy);
        EXPECT_EQ(dt.hour, 1);
        EXPECT_EQ(dt.minute, 1);
        EXPECT_EQ(dt.second, 1);
    }
}

// ---- the heap tripwire ----------------------------------------------------------------------------

TEST(Allocation, NoApiRoutineTouchesTheHeap) {
    // Everything in api.hpp is fixed-size arithmetic written through the caller's spans, so the
    // process-wide counter must not move across a run. Per alloc_counter.hpp the check with teeth
    // is the one around the SECOND run: a lane that allocated a workspace once and reused it would
    // slip past a check that only wrapped the first.
    const Rotations r = rotations_for(kYear, kDoy, kUt);
    const std::array<Rotations, 1> shared{r};
    std::vector<fx::vec3d> in(256, fx::vec3d{2.0, 1.0, -3.0});
    std::vector<fx::vec3d> out(256);
    std::vector<Status> st(256);
    std::vector<Position<Frame::GEO>> pts(256, Position<Frame::GEO>{fx::vec3d{2.0, 1.0, -3.0}});
    std::vector<double> mlt(256);
    std::vector<Position<Frame::GSM>> gsm(256);

    volatile double sink = 0.0;
    const auto exercise = [&]() {
        for (int a : kSysaxes) {
            const Result<bool> ct = api::coord_trans_vec(1, a, in, out, shared, st);
            sink = sink + out[0][0] + static_cast<double>(ct.value);
        }
        const Result<bool> m = api::get_mlt_vec(pts, shared, mlt, st);
        const Result<bool> t = api::transform_vec<Frame::GSM>(
            std::span<const Position<Frame::GEO>>{pts}, std::span<Position<Frame::GSM>>{gsm},
            shared);
        sink = sink + mlt[0] + gsm[0].v[0] + static_cast<double>(m.value) +
               static_cast<double>(t.value);
    };

    exercise();
    const std::size_t before = cheatah_space_test::allocation_count();
    exercise();
    EXPECT_EQ(cheatah_space_test::allocation_count(), before);
    EXPECT_TRUE(std::isfinite(static_cast<double>(sink)));
}
