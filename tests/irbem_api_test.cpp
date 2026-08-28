// Unit tests for space.irbem's IRBEM-named adapter layer — api.hpp.
//
// The layer has no physics of its own, so these tests ask three questions and nothing else:
//
//   1. **Does the name mean what IRBEM's name means?** Answered differentially, against golden
//      values harvested from the reference library as a BLACK BOX (dlopen + the documented C entry
//      points on /tmp/irbem-builds/libirbem-O2.so; see tools/oracle/). The literals below carry the
//      exact call that produced them, and are printed at %.17g so the literal round-trips to the
//      double the reference actually returned — at %.12g the truncation itself was 1e-12 relative,
//      which is above the residual of the legs that agree bit for bit. Nothing here links IRBEM, so
//      the QA gate can run this file on a machine that has never seen the reference.
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
//
// And on the SPATIAL sweep, which an earlier revision of this file did not have and needed. Eight
// epochs of one GEO point is eight rows of one geometry: it fixed a longitude (+165 deg), a radius
// (5.45 Re) and a quadrant, and each of those hid a real disagreement with the reference. The
// longitude hid that the reference folds SPH into [0, 360) while leaving GDZ and RLL in
// (-180, 180], so ->SPH read as EXACT when it was a full 360 deg out for every point with a
// negative GEO y. The radius hid that the iterative-versus-closed-form geodetic gap grows towards
// the surface, so the ->GDZ and ->RLL caps were set an order and a half tighter than the routine
// can hold. Neither is visible without a second point, so there are now four seeds spanning both
// signs of y, 1.04 to 5.45 Re and both hemispheres, plus a polar-axis test of its own — the axis
// being the one place the reference is WRONG and this module is not.
#include <gtest/gtest.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <string>
#include <iterator>
#include <limits>
#include <numbers>
#include <optional>
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
    {1965, 1, 0.0, fx::vec3d{6.6, 0, 0}, 23.630606468980492},
    {1965, 1, 0.0, fx::vec3d{0, 6.6, 0}, 5.580930424733018},
    {1965, 1, 0.0, fx::vec3d{-4.2, 1.1, 3.3}, 10.22334269395362},
    {1965, 1, 0.0, fx::vec3d{1.5, -1.5, 4}, 21.565814682180928},
    {1965, 1, 0.0, fx::vec3d{2, 3, -5}, 2.3199146839799667},
    {1980, 80, 3600.0, fx::vec3d{6.6, 0, 0}, 0.86142429022752154},
    {1980, 80, 3600.0, fx::vec3d{0, 6.6, 0}, 6.8159303231443777},
    {1980, 80, 3600.0, fx::vec3d{-4.2, 1.1, 3.3}, 11.459066258081901},
    {1980, 80, 3600.0, fx::vec3d{1.5, -1.5, 4}, 22.788283143732208},
    {1980, 80, 3600.0, fx::vec3d{2, 3, -5}, 3.5889248900161448},
    {1996, 241, 60360.0, fx::vec3d{6.6, 0, 0}, 16.767965034116401},
    {1996, 241, 60360.0, fx::vec3d{0, 6.6, 0}, 22.728246130364056},
    {1996, 241, 60360.0, fx::vec3d{-4.2, 1.1, 3.3}, 3.3816893051071304},
    {1996, 241, 60360.0, fx::vec3d{1.5, -1.5, 4}, 14.644901850419025},
    {1996, 241, 60360.0, fx::vec3d{2, 3, -5}, 19.562919964897805},
    {2003, 300, 70000.0, fx::vec3d{6.6, 0, 0}, 19.655573338725542},
    {2003, 300, 70000.0, fx::vec3d{0, 6.6, 0}, 1.6185358939650722},
    {2003, 300, 70000.0, fx::vec3d{-4.2, 1.1, 3.3}, 6.2783887110182608},
    {2003, 300, 70000.0, fx::vec3d{1.5, -1.5, 4}, 17.50417750320004},
    {2003, 300, 70000.0, fx::vec3d{2, 3, -5}, 22.484768398660155},
    {2015, 180, 43200.0, fx::vec3d{6.6, 0, 0}, 11.674330368781117},
    {2015, 180, 43200.0, fx::vec3d{0, 6.6, 0}, 17.643231682107746},
    {2015, 180, 43200.0, fx::vec3d{-4.2, 1.1, 3.3}, 22.316306992611864},
    {2015, 180, 43200.0, fx::vec3d{1.5, -1.5, 4}, 9.4686715980607818},
    {2015, 180, 43200.0, fx::vec3d{2, 3, -5}, 14.579943533910409},
    {2020, 90, 21600.0, fx::vec3d{6.6, 0, 0}, 5.9429690500656598},
    {2020, 90, 21600.0, fx::vec3d{0, 6.6, 0}, 11.913671666583237},
    {2020, 90, 21600.0, fx::vec3d{-4.2, 1.1, 3.3}, 16.594262388902031},
    {2020, 90, 21600.0, fx::vec3d{1.5, -1.5, 4}, 3.7084668677714596},
    {2020, 90, 21600.0, fx::vec3d{2, 3, -5}, 8.8768668438998475},
    {2024, 355, 86399.0, fx::vec3d{6.6, 0, 0}, 23.773069810563186},
    {2024, 355, 86399.0, fx::vec3d{0, 6.6, 0}, 5.7448747119229475},
    {2024, 355, 86399.0, fx::vec3d{-4.2, 1.1, 3.3}, 10.429919004721766},
    {2024, 355, 86399.0, fx::vec3d{1.5, -1.5, 4}, 21.521982999620967},
    {2024, 355, 86399.0, fx::vec3d{2, 3, -5}, 2.7241972470641933},
    {2029, 172, 7200.0, fx::vec3d{6.6, 0, 0}, 2.1316544718942261},
    {2029, 172, 7200.0, fx::vec3d{0, 6.6, 0}, 8.1049525973705823},
    {2029, 172, 7200.0, fx::vec3d{-4.2, 1.1, 3.3}, 12.795158808322975},
    {2029, 172, 7200.0, fx::vec3d{1.5, -1.5, 4}, 23.86198996911336},
    {2029, 172, 7200.0, fx::vec3d{2, 3, -5}, 5.1051944455405609},
};

// coord_trans goldens: coord_trans_vec1_(1,sysin,sysout,iyr,idoy,secs,xIN,xOUT)
constexpr CoordGolden kCoordGoldens[] = {
    {1965, 1, 0.0, 1, 0, fx::vec3d{-4.2, 1.1, 3.3}, fx::vec3d{28374.64362687297, 37.271646567743183, 165.32360686255001}},
    {1965, 1, 0.0, 1, 1, fx::vec3d{-4.2, 1.1, 3.3}, fx::vec3d{-4.2000000000000002, 1.1000000000000001, 3.2999999999999998}},
    {1965, 1, 0.0, 1, 2, fx::vec3d{-4.2, 1.1, 3.3}, fx::vec3d{2.5582103410202439, -2.1143772468053381, 4.3272356891309913}},
    {1965, 1, 0.0, 1, 3, fx::vec3d{-4.2, 1.1, 3.3}, fx::vec3d{2.5582103410202439, -0.7945955690106602, 4.7501366837052652}},
    {1965, 1, 0.0, 1, 4, fx::vec3d{-4.2, 1.1, 3.3}, fx::vec3d{4.2131523119428795, -2.1143772468053381, 2.742035093605157}},
    {1965, 1, 0.0, 1, 5, fx::vec3d{-4.2, 1.1, 3.3}, fx::vec3d{-0.32093738375299785, -4.3297805020242972, 3.2999999999999998}},
    {1965, 1, 0.0, 1, 6, fx::vec3d{-4.2, 1.1, 3.3}, fx::vec3d{-3.0835080747432309, -3.5655604746562717, 2.742035093605157}},
    {1965, 1, 0.0, 1, 7, fx::vec3d{-4.2, 1.1, 3.3}, fx::vec3d{5.4534392817743926, 37.237673077844335, 165.32360686255001}},
    {1965, 1, 0.0, 1, 8, fx::vec3d{-4.2, 1.1, 3.3}, fx::vec3d{5.4534392817743926, 37.271646567743183, 165.32360686255001}},
    {1980, 80, 3600.0, 1, 0, fx::vec3d{-4.2, 1.1, 3.3}, fx::vec3d{28374.64362687297, 37.271646567743183, 165.32360686255001}},
    {1980, 80, 3600.0, 1, 1, fx::vec3d{-4.2, 1.1, 3.3}, fx::vec3d{-4.2000000000000002, 1.1000000000000001, 3.2999999999999998}},
    {1980, 80, 3600.0, 1, 2, fx::vec3d{-4.2, 1.1, 3.3}, fx::vec3d{4.3304264894157187, -0.66320086714467186, 3.2477024231886169}},
    {1980, 80, 3600.0, 1, 3, fx::vec3d{-4.2, 1.1, 3.3}, fx::vec3d{4.3304264894157187, 1.2084152611244356, 3.0866963408422299}},
    {1980, 80, 3600.0, 1, 4, fx::vec3d{-4.2, 1.1, 3.3}, fx::vec3d{4.6517401475689235, -0.66320086714467186, 2.7679375371048418}},
    {1980, 80, 3600.0, 1, 5, fx::vec3d{-4.2, 1.1, 3.3}, fx::vec3d{4.3391292256226563, -0.14818084676277066, 3.2999999999999998}},
    {1980, 80, 3600.0, 1, 6, fx::vec3d{-4.2, 1.1, 3.3}, fx::vec3d{-3.015548846696448, -3.6034687377405827, 2.7679375371048418}},
    {1980, 80, 3600.0, 1, 7, fx::vec3d{-4.2, 1.1, 3.3}, fx::vec3d{5.4534392817743926, 37.237673077844335, 165.32360686255001}},
    {1980, 80, 3600.0, 1, 8, fx::vec3d{-4.2, 1.1, 3.3}, fx::vec3d{5.4534392817743926, 37.271646567743183, 165.32360686255001}},
    {1996, 241, 60360.0, 1, 0, fx::vec3d{-4.2, 1.1, 3.3}, fx::vec3d{28374.64362687297, 37.271646567743183, 165.32360686255001}},
    {1996, 241, 60360.0, 1, 1, fx::vec3d{-4.2, 1.1, 3.3}, fx::vec3d{-4.2000000000000002, 1.1000000000000001, 3.2999999999999998}},
    {1996, 241, 60360.0, 1, 2, fx::vec3d{-4.2, 1.1, 3.3}, fx::vec3d{-1.8188513069118861, -3.6200676868204433, 3.6506013020026451}},
    {1996, 241, 60360.0, 1, 3, fx::vec3d{-4.2, 1.1, 3.3}, fx::vec3d{-1.8188513069118861, -4.7097627964717743, 2.061864383831701}},
    {1996, 241, 60360.0, 1, 4, fx::vec3d{-4.2, 1.1, 3.3}, fx::vec3d{-2.9603184710242423, -3.6200676868204433, 2.8056415474809278}},
    {1996, 241, 60360.0, 1, 5, fx::vec3d{-4.2, 1.1, 3.3}, fx::vec3d{3.5992337745362475, 2.4280684167126667, 3.2999999999999998}},
    {1996, 241, 60360.0, 1, 6, fx::vec3d{-4.2, 1.1, 3.3}, fx::vec3d{-2.9452604144816696, -3.6323293625352986, 2.8056415474809278}},
    {1996, 241, 60360.0, 1, 7, fx::vec3d{-4.2, 1.1, 3.3}, fx::vec3d{5.4534392817743926, 37.237673077844335, 165.32360686255001}},
    {1996, 241, 60360.0, 1, 8, fx::vec3d{-4.2, 1.1, 3.3}, fx::vec3d{5.4534392817743926, 37.271646567743183, 165.32360686255001}},
    {2003, 300, 70000.0, 1, 0, fx::vec3d{-4.2, 1.1, 3.3}, fx::vec3d{28374.64362687297, 37.271646567743183, 165.32360686255001}},
    {2003, 300, 70000.0, 1, 1, fx::vec3d{-4.2, 1.1, 3.3}, fx::vec3d{-4.2000000000000002, 1.1000000000000001, 3.2999999999999998}},
    {2003, 300, 70000.0, 1, 2, fx::vec3d{-4.2, 1.1, 3.3}, fx::vec3d{0.075660497065418775, -4.6529200224943601, 2.8433449937450259}},
    {2003, 300, 70000.0, 1, 3, fx::vec3d{-4.2, 1.1, 3.3}, fx::vec3d{0.075660497065418775, -5.1615054191481518, 1.7591645170552961}},
    {2003, 300, 70000.0, 1, 4, fx::vec3d{-4.2, 1.1, 3.3}, fx::vec3d{0.33971580191380929, -4.6529200224943601, 2.8239915789890633}},
    {2003, 300, 70000.0, 1, 5, fx::vec3d{-4.2, 1.1, 3.3}, fx::vec3d{-2.9462278306850314, 3.1890032250998082, 3.2999999999999998}},
    {2003, 300, 70000.0, 1, 6, fx::vec3d{-4.2, 1.1, 3.3}, fx::vec3d{-2.913660540524929, -3.6435770633248361, 2.8239915789890633}},
    {2003, 300, 70000.0, 1, 7, fx::vec3d{-4.2, 1.1, 3.3}, fx::vec3d{5.4534392817743926, 37.237673077844335, 165.32360686255001}},
    {2003, 300, 70000.0, 1, 8, fx::vec3d{-4.2, 1.1, 3.3}, fx::vec3d{5.4534392817743926, 37.271646567743183, 165.32360686255001}},
    {2015, 180, 43200.0, 1, 0, fx::vec3d{-4.2, 1.1, 3.3}, fx::vec3d{28374.64362687297, 37.271646567743183, 165.32360686255001}},
    {2015, 180, 43200.0, 1, 1, fx::vec3d{-4.2, 1.1, 3.3}, fx::vec3d{-4.2000000000000002, 1.1000000000000001, 3.2999999999999998}},
    {2015, 180, 43200.0, 1, 2, fx::vec3d{-4.2, 1.1, 3.3}, fx::vec3d{-2.5425554394306351, 1.9793657824267803, 4.3997185065387665}},
    {2015, 180, 43200.0, 1, 3, fx::vec3d{-4.2, 1.1, 3.3}, fx::vec3d{-2.5425554394306351, 0.89807325698354279, 4.7403841987319062}},
    {2015, 180, 43200.0, 1, 4, fx::vec3d{-4.2, 1.1, 3.3}, fx::vec3d{-4.1958318871635205, 1.9793657824267803, 2.8665494717551674}},
    {2015, 180, 43200.0, 1, 5, fx::vec3d{-4.2, 1.1, 3.3}, fx::vec3d{-0.5609795043219894, -4.3052644513119818, 3.2999999999999998}},
    {2015, 180, 43200.0, 1, 6, fx::vec3d{-4.2, 1.1, 3.3}, fx::vec3d{-2.8254425657542344, -3.6796424328464719, 2.8665494717551674}},
    {2015, 180, 43200.0, 1, 7, fx::vec3d{-4.2, 1.1, 3.3}, fx::vec3d{5.4534392817743926, 37.237673077844335, 165.32360686255001}},
    {2015, 180, 43200.0, 1, 8, fx::vec3d{-4.2, 1.1, 3.3}, fx::vec3d{5.4534392817743926, 37.271646567743183, 165.32360686255001}},
    {2020, 90, 21600.0, 1, 0, fx::vec3d{-4.2, 1.1, 3.3}, fx::vec3d{28374.64362687297, 37.271646567743183, 165.32360686255001}},
    {2020, 90, 21600.0, 1, 1, fx::vec3d{-4.2, 1.1, 3.3}, fx::vec3d{-4.2000000000000002, 1.1000000000000001, 3.2999999999999998}},
    {2020, 90, 21600.0, 1, 2, fx::vec3d{-4.2, 1.1, 3.3}, fx::vec3d{1.4054550334981459, 4.3206359065262143, 3.01609043499215}},
    {2020, 90, 21600.0, 1, 3, fx::vec3d{-4.2, 1.1, 3.3}, fx::vec3d{1.4054550334981459, 5.1033943596596503, 1.3120374261584546}},
    {2020, 90, 21600.0, 1, 4, fx::vec3d{-4.2, 1.1, 3.3}, fx::vec3d{1.6659873929798863, 4.3206359065262143, 2.8803804209979784}},
    {2020, 90, 21600.0, 1, 5, fx::vec3d{-4.2, 1.1, 3.3}, fx::vec3d{0.49800155230881293, 4.3130029508334466, 3.2999999999999998}},
    {2020, 90, 21600.0, 1, 6, fx::vec3d{-4.2, 1.1, 3.3}, fx::vec3d{-2.8078440503943218, -3.6823118307656562, 2.8803804209979784}},
    {2020, 90, 21600.0, 1, 7, fx::vec3d{-4.2, 1.1, 3.3}, fx::vec3d{5.4534392817743926, 37.237673077844335, 165.32360686255001}},
    {2020, 90, 21600.0, 1, 8, fx::vec3d{-4.2, 1.1, 3.3}, fx::vec3d{5.4534392817743926, 37.271646567743183, 165.32360686255001}},
    {2024, 355, 86399.0, 1, 0, fx::vec3d{-4.2, 1.1, 3.3}, fx::vec3d{28374.64362687297, 37.271646567743183, 165.32360686255001}},
    {2024, 355, 86399.0, 1, 1, fx::vec3d{-4.2, 1.1, 3.3}, fx::vec3d{-4.2000000000000002, 1.1000000000000001, 3.2999999999999998}},
    {2024, 355, 86399.0, 1, 2, fx::vec3d{-4.2, 1.1, 3.3}, fx::vec3d{2.5494053747536869, -1.8481373265080876, 4.4525184623474114}},
    {2024, 355, 86399.0, 1, 3, fx::vec3d{-4.2, 1.1, 3.3}, fx::vec3d{2.5494053747536869, -1.0787377684342689, 4.6989083718049001}},
    {2024, 355, 86399.0, 1, 4, fx::vec3d{-4.2, 1.1, 3.3}, fx::vec3d{4.2400563300442338, -1.8481373265080876, 2.888998224371099}},
    {2024, 355, 86399.0, 1, 5, fx::vec3d{-4.2, 1.1, 3.3}, fx::vec3d{-1.0961692035500177, -4.20100143741805, 3.2999999999999998}},
    {2024, 355, 86399.0, 1, 6, fx::vec3d{-4.2, 1.1, 3.3}, fx::vec3d{-2.7953675774083262, -3.6850521525149333, 2.888998224371099}},
    {2024, 355, 86399.0, 1, 7, fx::vec3d{-4.2, 1.1, 3.3}, fx::vec3d{5.4534392817743926, 37.237673077844335, 165.32360686255001}},
    {2024, 355, 86399.0, 1, 8, fx::vec3d{-4.2, 1.1, 3.3}, fx::vec3d{5.4534392817743926, 37.271646567743183, 165.32360686255001}},
    {2029, 172, 7200.0, 1, 0, fx::vec3d{-4.2, 1.1, 3.3}, fx::vec3d{28374.64362687297, 37.271646567743183, 165.32360686255001}},
    {2029, 172, 7200.0, 1, 1, fx::vec3d{-4.2, 1.1, 3.3}, fx::vec3d{-4.2000000000000002, 1.1000000000000001, 3.2999999999999998}},
    {2029, 172, 7200.0, 1, 2, fx::vec3d{-4.2, 1.1, 3.3}, fx::vec3d{5.1625078089833698, 0.95439573992793902, 1.4756835344317996}},
    {2029, 172, 7200.0, 1, 3, fx::vec3d{-4.2, 1.1, 3.3}, fx::vec3d{5.1625078089833698, 1.114381350651229, 1.3585613188217756}},
    {2029, 172, 7200.0, 1, 4, fx::vec3d{-4.2, 1.1, 3.3}, fx::vec3d{4.5182296312527805, 0.95439573992793902, 2.9008153631137512}},
    {2029, 172, 7200.0, 1, 5, fx::vec3d{-4.2, 1.1, 3.3}, fx::vec3d{-1.1149729071177692, 4.1960499778235905, 3.2999999999999998}},
    {2029, 172, 7200.0, 1, 6, fx::vec3d{-4.2, 1.1, 3.3}, fx::vec3d{-2.7732911772063153, -3.6924417768670157, 2.9008153631137512}},
    {2029, 172, 7200.0, 1, 7, fx::vec3d{-4.2, 1.1, 3.3}, fx::vec3d{5.4534392817743926, 37.237673077844335, 165.32360686255001}},
    {2029, 172, 7200.0, 1, 8, fx::vec3d{-4.2, 1.1, 3.3}, fx::vec3d{5.4534392817743926, 37.271646567743183, 165.32360686255001}},
    {2015, 180, 43200.0, 0, 1, fx::vec3d{500, 33, -75}, fx::vec3d{0.23455110640987956, -0.87535664609316477, 0.58486284506188879}},
    {2015, 180, 43200.0, 0, 2, fx::vec3d{500, 33, -75}, fx::vec3d{0.43417178828460656, -0.78434588742437694, 0.59968693003072293}},
    {2015, 180, 43200.0, 0, 4, fx::vec3d{500, 33, -75}, fx::vec3d{0.13197252881770383, -0.78434588742437694, 0.72850038258698713}},
    {2015, 180, 43200.0, 7, 1, fx::vec3d{1.3, 33, -75}, fx::vec3d{0.28218309021659038, -1.0531216297251076, 0.70803074551953538}},
    {2015, 180, 43200.0, 7, 2, fx::vec3d{1.3, 33, -75}, fx::vec3d{0.52407575828112662, -0.94290793816109253, 0.72544415342248048}},
    {2015, 180, 43200.0, 7, 4, fx::vec3d{1.3, 33, -75}, fx::vec3d{0.15861647630055134, -0.94290793816109253, 0.88077547286399493}},
    {2015, 180, 43200.0, 8, 1, fx::vec3d{1.3, 33, -75}, fx::vec3d{0.28261424270647728, -1.0547307127231753, 0.70545893817989658}},
    {2015, 180, 43200.0, 8, 2, fx::vec3d{1.3, 33, -75}, fx::vec3d{0.52343555168260081, -0.94494768929774831, 0.72324898045249775}},

    // Three MORE GEO seeds, at 1965, 2015 and 2029. The eight-epoch block above is one spatial
    // point, {-4.2, 1.1, 3.3}, and a single point cannot see a convention: its GEO longitude is
    // +165 deg, so it never exercises the [0, 360) fold the reference applies to SPH and not to
    // GDZ or RLL, and it sits at 5.45 Re, where the iterative-versus-closed-form geodetic gap is
    // two orders below what it is near the surface. {1.5, -1.5, 4} is the negative-y quadrant
    // (SPH longitude 315 against RLL's -45), {0.62, -0.38, 0.75} is 1.04 Re, and
    // {-0.9, -0.7, -1.05} is southern, negative-y AND near the surface at 1.55 Re.
    {1965, 1, 0, 1, 0, fx::vec3d{1.5, -1.5, 4}, fx::vec3d{22485.387632972968, 62.096812811799659, -45}},
    {1965, 1, 0, 1, 1, fx::vec3d{1.5, -1.5, 4}, fx::vec3d{1.5, -1.5, 4}},
    {1965, 1, 0, 1, 2, fx::vec3d{1.5, -1.5, 4}, fx::vec3d{-2.925235329750862, 0.8388787563490695, 3.3525036461909488}},
    {1965, 1, 0, 1, 3, fx::vec3d{1.5, -1.5, 4}, fx::vec3d{-2.925235329750862, 1.758877541736096, 2.9750761914529402}},
    {1965, 1, 0, 1, 4, fx::vec3d{1.5, -1.5, 4}, fx::vec3d{-1.1331508933769661, 0.8388787563490695, 4.3025866040075496}},
    {1965, 1, 0, 1, 5, fx::vec3d{1.5, -1.5, 4}, fx::vec3d{1.2034446437798045, 1.7469175679916036, 4}},
    {1965, 1, 0, 1, 6, fx::vec3d{1.5, -1.5, 4}, fx::vec3d{1.0912718483325228, 0.89267814359470987, 4.3025866040075504}},
    {1965, 1, 0, 1, 7, fx::vec3d{1.5, -1.5, 4}, fx::vec3d{4.5276925690687087, 62.06164727039765, 315}},
    {1965, 1, 0, 1, 8, fx::vec3d{1.5, -1.5, 4}, fx::vec3d{4.5276925690687087, 62.096812811799659, -45}},
    {1965, 1, 0, 1, 0, fx::vec3d{0.62, -0.38, 0.75}, fx::vec3d{288.60440761699033, 46.068779058089078, -31.504266719204189}},
    {1965, 1, 0, 1, 1, fx::vec3d{0.62, -0.38, 0.75}, fx::vec3d{0.62, -0.38, 0.75}},
    {1965, 1, 0, 1, 2, fx::vec3d{0.62, -0.38, 0.75}, fx::vec3d{-0.85883877746938708, 0.2863547694908547, 0.52124552785105815}},
    {1965, 1, 0, 1, 3, fx::vec3d{0.62, -0.38, 0.75}, fx::vec3d{-0.85883877746938708, 0.42294361041218836, 0.41818530814066485}},
    {1965, 1, 0, 1, 4, fx::vec3d{0.62, -0.38, 0.75}, fx::vec3d{-0.53822257044698341, 0.2863547694908547, 0.84830266453152292}},
    {1965, 1, 0, 1, 5, fx::vec3d{0.62, -0.38, 0.75}, fx::vec3d{0.26139480915393992, 0.67858142749958561, 0.75}},
    {1965, 1, 0, 1, 6, fx::vec3d{0.62, -0.38, 0.75}, fx::vec3d{0.40966281419980138, 0.45150743959609085, 0.84830266453152292}},
    {1965, 1, 0, 1, 7, fx::vec3d{0.62, -0.38, 0.75}, fx::vec3d{1.0446530524533013, 45.884804721120098, 328.49573328079583}},
    {1965, 1, 0, 1, 8, fx::vec3d{0.62, -0.38, 0.75}, fx::vec3d{1.0446530524533013, 46.068779058089078, -31.504266719204189}},
    {1965, 1, 0, 1, 0, fx::vec3d{-0.90000000000000002, -0.69999999999999996, -1.05}, fx::vec3d{3507.0503343306423, -42.765987670271478, -142.1250163489018}},
    {1965, 1, 0, 1, 1, fx::vec3d{-0.90000000000000002, -0.69999999999999996, -1.05}, fx::vec3d{-0.90000000000000002, -0.69999999999999996, -1.05}},
    {1965, 1, 0, 1, 2, fx::vec3d{-0.90000000000000002, -0.69999999999999996, -1.05}, fx::vec3d{1.2485503913153186, 0.79911369533051402, -0.45281256859938834}},
    {1965, 1, 0, 1, 3, fx::vec3d{-0.90000000000000002, -0.69999999999999996, -1.05}, fx::vec3d{1.2485503913153186, 0.63709716433875196, -0.6617239304725393}},
    {1965, 1, 0, 1, 4, fx::vec3d{-0.90000000000000002, -0.69999999999999996, -1.05}, fx::vec3d{0.91781373077385331, 0.79911369533051402, -0.95996627937557844}},
    {1965, 1, 0, 1, 5, fx::vec3d{-0.90000000000000002, -0.69999999999999996, -1.05}, fx::vec3d{0.8514597266768682, -0.75829831454866936, -1.05}},
    {1965, 1, 0, 1, 6, fx::vec3d{-0.90000000000000002, -0.69999999999999996, -1.05}, fx::vec3d{0.54948433420209797, -1.0858322655586756, -0.95996627937557844}},
    {1965, 1, 0, 1, 7, fx::vec3d{-0.90000000000000002, -0.69999999999999996, -1.05}, fx::vec3d{1.55, -42.642309979939938, 217.8749836510982}},
    {1965, 1, 0, 1, 8, fx::vec3d{-0.90000000000000002, -0.69999999999999996, -1.05}, fx::vec3d{1.55, -42.765987670271478, -142.1250163489018}},
    {2015, 180, 43200, 1, 0, fx::vec3d{1.5, -1.5, 4}, fx::vec3d{22485.387632972968, 62.096812811799659, -45}},
    {2015, 180, 43200, 1, 1, fx::vec3d{1.5, -1.5, 4}, fx::vec3d{1.5, -1.5, 4}},
    {2015, 180, 43200, 1, 2, fx::vec3d{1.5, -1.5, 4}, fx::vec3d{2.9352454330823274, -0.94590442070122449, 3.3150564210081286}},
    {2015, 180, 43200, 1, 3, fx::vec3d{1.5, -1.5, 4}, fx::vec3d{2.9352454330823274, -1.6934000630214483, 3.0027448900684899}},
    {2015, 180, 43200, 1, 4, fx::vec3d{1.5, -1.5, 4}, fx::vec3d{1.2120132086448443, -0.94590442070122449, 4.2586721884841419}},
    {2015, 180, 43200, 1, 5, fx::vec3d{1.5, -1.5, 4}, fx::vec3d{1.2986348412711155, 1.6773632728293368, 4}},
    {2015, 180, 43200, 1, 6, fx::vec3d{1.5, -1.5, 4}, fx::vec3d{1.1817653465526656, 0.98343370733311175, 4.2586721884841419}},
    {2015, 180, 43200, 1, 7, fx::vec3d{1.5, -1.5, 4}, fx::vec3d{4.5276925690687087, 62.06164727039765, 315}},
    {2015, 180, 43200, 1, 8, fx::vec3d{1.5, -1.5, 4}, fx::vec3d{4.5276925690687087, 62.096812811799659, -45}},
    {2015, 180, 43200, 1, 0, fx::vec3d{0.62, -0.38, 0.75}, fx::vec3d{288.60440761699033, 46.068779058089078, -31.504266719204189}},
    {2015, 180, 43200, 1, 1, fx::vec3d{0.62, -0.38, 0.75}, fx::vec3d{0.62, -0.38, 0.75}},
    {2015, 180, 43200, 1, 2, fx::vec3d{0.62, -0.38, 0.75}, fx::vec3d{0.86026126484819154, -0.30321850892090463, 0.50922400969472448}},
    {2015, 180, 43200, 1, 3, fx::vec3d{0.62, -0.38, 0.75}, fx::vec3d{0.86026126484819154, -0.41368516724908844, 0.42437887483421793}},
    {2015, 180, 43200, 1, 4, fx::vec3d{0.62, -0.38, 0.75}, fx::vec3d{0.55525657414200602, -0.30321850892090463, 0.83129337343675869}},
    {2015, 180, 43200, 1, 5, fx::vec3d{0.62, -0.38, 0.75}, fx::vec3d{0.29868921859735831, 0.66301187824480157, 0.75}},
    {2015, 180, 43200, 1, 6, fx::vec3d{0.62, -0.38, 0.75}, fx::vec3d{0.4142423439044246, 0.47817842673703098, 0.83129337343675869}},
    {2015, 180, 43200, 1, 7, fx::vec3d{0.62, -0.38, 0.75}, fx::vec3d{1.0446530524533013, 45.884804721120098, 328.49573328079583}},
    {2015, 180, 43200, 1, 8, fx::vec3d{0.62, -0.38, 0.75}, fx::vec3d{1.0446530524533013, 46.068779058089078, -31.504266719204189}},
    {2015, 180, 43200, 1, 0, fx::vec3d{-0.90000000000000002, -0.69999999999999996, -1.05}, fx::vec3d{3507.0503343306423, -42.765987670271478, -142.1250163489018}},
    {2015, 180, 43200, 1, 1, fx::vec3d{-0.90000000000000002, -0.69999999999999996, -1.05}, fx::vec3d{-0.90000000000000002, -0.69999999999999996, -1.05}},
    {2015, 180, 43200, 1, 2, fx::vec3d{-0.90000000000000002, -0.69999999999999996, -1.05}, fx::vec3d{-1.250668186693207, -0.78357343643047028, -0.47364729125591881}},
    {2015, 180, 43200, 1, 3, fx::vec3d{-0.90000000000000002, -0.69999999999999996, -1.05}, fx::vec3d{-1.250668186693207, -0.65143270525446084, -0.64337295808048611}},
    {2015, 180, 43200, 1, 4, fx::vec3d{-0.90000000000000002, -0.69999999999999996, -1.05}, fx::vec3d{-0.92262380406225675, -0.78357343643047028, -0.96813118217431382}},
    {2015, 180, 43200, 1, 5, fx::vec3d{-0.90000000000000002, -0.69999999999999996, -1.05}, fx::vec3d{0.80801808942423858, -0.80442946686655081, -1.05}},
    {2015, 180, 43200, 1, 6, fx::vec3d{-0.90000000000000002, -0.69999999999999996, -1.05}, fx::vec3d{0.56971594221413602, -1.0680101868843879, -0.96813118217431382}},
    {2015, 180, 43200, 1, 7, fx::vec3d{-0.90000000000000002, -0.69999999999999996, -1.05}, fx::vec3d{1.55, -42.642309979939938, 217.8749836510982}},
    {2015, 180, 43200, 1, 8, fx::vec3d{-0.90000000000000002, -0.69999999999999996, -1.05}, fx::vec3d{1.55, -42.765987670271478, -142.1250163489018}},
    {2029, 172, 7200, 1, 0, fx::vec3d{1.5, -1.5, 4}, fx::vec3d{22485.387632972968, 62.096812811799659, -45}},
    {2029, 172, 7200, 1, 1, fx::vec3d{1.5, -1.5, 4}, fx::vec3d{1.5, -1.5, 4}},
    {2029, 172, 7200, 1, 2, fx::vec3d{1.5, -1.5, 4}, fx::vec3d{-0.28513071213305863, 0.056935913564766572, 4.518346907746797}},
    {2029, 172, 7200, 1, 3, fx::vec3d{1.5, -1.5, 4}, fx::vec3d{-0.28513071213305863, 0.56490793157897101, 4.4833948779716772}},
    {2029, 172, 7200, 1, 4, fx::vec3d{1.5, -1.5, 4}, fx::vec3d{-1.5751356960872978, 0.056935913564766572, 4.2444912346072909}},
    {2029, 172, 7200, 1, 5, fx::vec3d{1.5, -1.5, 4}, fx::vec3d{-0.56485929895455855, -2.0447332276814412, 4}},
    {2029, 172, 7200, 1, 6, fx::vec3d{1.5, -1.5, 4}, fx::vec3d{1.2232935429154415, 0.99390495883824759, 4.2444912346072918}},
    {2029, 172, 7200, 1, 7, fx::vec3d{1.5, -1.5, 4}, fx::vec3d{4.5276925690687087, 62.06164727039765, 315}},
    {2029, 172, 7200, 1, 8, fx::vec3d{1.5, -1.5, 4}, fx::vec3d{4.5276925690687087, 62.096812811799659, -45}},
    {2029, 172, 7200, 1, 0, fx::vec3d{0.62, -0.38, 0.75}, fx::vec3d{288.60440761699033, 46.068779058089078, -31.504266719204189}},
    {2029, 172, 7200, 1, 1, fx::vec3d{0.62, -0.38, 0.75}, fx::vec3d{0.62, -0.38, 0.75}},
    {2029, 172, 7200, 1, 2, fx::vec3d{0.62, -0.38, 0.75}, fx::vec3d{-0.3685250107296687, -0.085338709347011238, 0.97375901595501746}},
    {2029, 172, 7200, 1, 3, fx::vec3d{0.62, -0.38, 0.75}, fx::vec3d{-0.3685250107296687, 0.024752657229611719, 0.97723270832048892}},
    {2029, 172, 7200, 1, 4, fx::vec3d{0.62, -0.38, 0.75}, fx::vec3d{-0.63350946857343149, -0.085338709347011238, 0.82624636635497217}},
    {2029, 172, 7200, 1, 5, fx::vec3d{0.62, -0.38, 0.75}, fx::vec3d{-0.024707774770337532, -0.72676648647684505, 0.75}},
    {2029, 172, 7200, 1, 6, fx::vec3d{0.62, -0.38, 0.75}, fx::vec3d{0.42075761888005192, 0.48122756388178145, 0.82624636635497217}},
    {2029, 172, 7200, 1, 7, fx::vec3d{0.62, -0.38, 0.75}, fx::vec3d{1.0446530524533013, 45.884804721120098, 328.49573328079583}},
    {2029, 172, 7200, 1, 8, fx::vec3d{0.62, -0.38, 0.75}, fx::vec3d{1.0446530524533013, 46.068779058089078, -31.504266719204189}},
    {2029, 172, 7200, 1, 0, fx::vec3d{-0.90000000000000002, -0.69999999999999996, -1.05}, fx::vec3d{3507.0503343306423, -42.765987670271478, -142.1250163489018}},
    {2029, 172, 7200, 1, 1, fx::vec3d{-0.90000000000000002, -0.69999999999999996, -1.05}, fx::vec3d{-0.90000000000000002, -0.69999999999999996, -1.05}},
    {2029, 172, 7200, 1, 2, fx::vec3d{-0.90000000000000002, -0.69999999999999996, -1.05}, fx::vec3d{-0.015997020352793823, 1.1741453572616423, -1.0117444219568317}},
    {2029, 172, 7200, 1, 3, fx::vec3d{-0.90000000000000002, -0.69999999999999996, -1.05}, fx::vec3d{-0.015997020352793823, 1.0528959331297858, -1.1374418740643111}},
    {2029, 172, 7200, 1, 4, fx::vec3d{-0.90000000000000002, -0.69999999999999996, -1.05}, fx::vec3d{0.27624707731817449, 1.1741453572616423, -0.97343219193434194}},
    {2029, 172, 7200, 1, 5, fx::vec3d{-0.90000000000000002, -0.69999999999999996, -1.05}, fx::vec3d{-1.0528671014997981, 0.4375738412878607, -1.05}},
    {2029, 172, 7200, 1, 6, fx::vec3d{-0.90000000000000002, -0.69999999999999996, -1.05}, fx::vec3d{0.56486435618618236, -1.0657664034938765, -0.97343219193434194}},
    {2029, 172, 7200, 1, 7, fx::vec3d{-0.90000000000000002, -0.69999999999999996, -1.05}, fx::vec3d{1.55, -42.642309979939938, 217.8749836510982}},
    {2029, 172, 7200, 1, 8, fx::vec3d{-0.90000000000000002, -0.69999999999999996, -1.05}, fx::vec3d{1.55, -42.765987670271478, -142.1250163489018}},
};

/// One `make_lstar1_(1, kext=0, options={1,0,0,0,0}, sysaxes=1, ...)` observation: the six outputs
/// IRBEM returns for one point, at its DEFAULT drift-shell resolution.
struct LstarGolden {
    int year;                        ///< `iyear`.
    int doy;                         ///< `idoy`.
    double ut;                       ///< `UT`.
    fx::vec3d geo;                   ///< the point, GEO, Earth radii.
    api::MagneticCoordinates out;    ///< `Lm`, `Lstar`, `Blocal`, `Bmin`, `XJ`, `MLT`.
};

constexpr LstarGolden kLstarGoldens[] = {
    {2015, 180, 43200.0, fx::vec3d{4, 0, 0}, {4.0633475080066876, 4.0484118075563575, 445.41161079063284, 444.85277462512221, 0.0021390145357650697, 11.674330368781117}},
    {2015, 180, 43200.0, fx::vec3d{0, 5, 0}, {5.0467734246650213, 5.031501108998266, 256.85346474122531, 231.88299758439467, 0.37201299706406205, 17.643231682107746}},
    {2015, 180, 43200.0, fx::vec3d{-6.6, 0, 0}, {6.5658181815039214, 6.5285816200395193, 107.39561820893306, 105.39905167993251, 0.086605938775183383, 23.674330368781114}},
    {2015, 180, 43200.0, fx::vec3d{3, 2, 1}, {3.9040976422198117, 3.8856478641774506, 604.45535880748002, 500.7768455815982, 0.52906929417203141, 14.058962009316195}},
    {2015, 180, 43200.0, fx::vec3d{-2, 4, -3}, {10.247147107956714, -9.9999999999999996e+30, 307.54745256517435, 27.661740220704633, 14.253752698425423, 19.518033518831757}},
    {2003, 300, 70000.0, fx::vec3d{4.5, -1, 0.5}, {4.8401945848776231, 4.8381319896781934, 301.41948320373859, 264.74704537601554, 0.45503169652254744, 18.894066399005762}},
    {2020, 90, 21600.0, fx::vec3d{-3, -3, 2}, {6.1081022807588283, 6.0724994551703535, 380.81031658558709, 130.40719752945904, 4.3631763800944183, 20.636381139255533}},
    {1996, 241, 60360.0, fx::vec3d{2.5, 1, -1}, {3.4261884893730947, 3.4364797204141198, 1424.3670794879015, 749.20009673009326, 1.531221287922667, 17.936540277322784}},
};


// ---- helpers -----------------------------------------------------------------------------------

/// The IGRF model for a (year, doy, ut) epoch — the argument @ref api::rotations_at needs.
/// An epoch outside IGRF-14's window is a fixture bug rather than a result, so it aborts the test.
[[nodiscard]] std::optional<ib::Igrf<>> model_at(int year, int doy, double ut) {
    const ib::DateTime dt = api::doy_and_ut2date_and_time(year, doy, ut);
    return ib::Igrf<>::at(
        api::date_and_time2decy(dt.year, dt.month, dt.day, dt.hour, dt.minute, dt.second));
}

/// The rotations for a (year, doy, ut) epoch with the IGRF sampled at the MIDDLE of the year —
/// which is the reference library's own convention, established by measurement rather than assumed.
///
/// GEO→MAG is a pure dipole rotation: no Sun, no sidereal time, only the degree-1 IGRF
/// coefficients. Sweeping day-of-year through 2015 and comparing GEO→MAG against
/// `coord_trans_vec1_` gives a residual that is 1.8e-3 Re at doy 1, falls to 3e-5 at doy 180 and
/// climbs symmetrically back to 1.7e-3 at doy 364 — a straight line through zero at mid-year. Pin
/// this module's IGRF to `year + 0.5` and the same comparison collapses to 1.3e-15, i.e. bit level.
/// So the reference samples the secular variation once per year, at the year's midpoint, and a
/// differential test that feeds it the exact-date coefficients is measuring that convention rather
/// than either library's arithmetic. The goldens below are therefore compared at MATCHED options.
[[nodiscard]] Rotations rotations_matched(int year, int doy, double ut) {
    const std::optional<ib::Igrf<>> m = ib::Igrf<>::at(static_cast<double>(year) + 0.5);
    EXPECT_TRUE(m.has_value()) << year;
    if (!m.has_value()) return Rotations{};
    const Result<Rotations> r = api::rotations_at(year, doy, ut, *m);
    EXPECT_EQ(r.status, Status::Ok);
    return r.value;
}

/// The rotations for a (year, doy, ut) epoch, built through the adapter's own entry point.
[[nodiscard]] Rotations rotations_for(int year, int doy, double ut) {
    const std::optional<ib::Igrf<>> m = model_at(year, doy, ut);
    EXPECT_TRUE(m.has_value()) << year;
    if (!m.has_value()) return Rotations{};
    const Result<Rotations> r = api::rotations_at(year, doy, ut, *m);
    EXPECT_EQ(r.status, Status::Ok);
    return r.value;
}

/// Pins the device lane OFF for the life of the object, and puts the environment back.
///
/// Needed because which lane a drift-shell batch runs on depends on its SIZE, so on a machine with
/// a device `make_lstar_vec` over n points and n separate `make_lstar` calls are not the same
/// arithmetic: n = 1 is 25 field lines and stays on the host, n = 6 is 150, which crosses
/// `irbem_igrf_f32`'s measured crossover of 128 and puts the flux quadrature in fp32 on the device.
/// Measured, that moves L* by ~1e-11 — far below every physics budget and far above `==`. So every
/// test that compares the batch lane to the scalar lane EXACTLY holds both to the host lane, and
/// the device-versus-host question is asked separately, with the tolerance it deserves, by
/// IrbemApi.MakeLstarVecDeviceLaneAgreesWithTheHostLane.
///
/// `CHEATAH_SPACE_IRBEM_NO_GPU` is the knob gpu/dispatch.hpp documents for exactly this, and it is
/// read fresh on every dispatch, so setting it mid-process works.
class HostLanePin {
public:
    HostLanePin() {
        if (const char* prev = std::getenv(kVar)) {
            had_ = true;
            prev_ = prev;
        }
        ::setenv(kVar, "1", 1);
    }
    HostLanePin(const HostLanePin&) = delete;
    HostLanePin& operator=(const HostLanePin&) = delete;
    HostLanePin(HostLanePin&&) = delete;
    HostLanePin& operator=(HostLanePin&&) = delete;
    ~HostLanePin() {
        if (had_) {
            ::setenv(kVar, prev_.c_str(), 1);
        } else {
            ::unsetenv(kVar);
        }
    }

private:
    static constexpr const char* kVar = "CHEATAH_SPACE_IRBEM_NO_GPU";
    bool had_ = false;
    std::string prev_;
};

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
    const std::optional<ib::Igrf<>> model = model_at(kYear, kDoy, kUt);
    ASSERT_TRUE(model.has_value());
    const ib::Igrf<>& m = *model;
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

    // An epoch IGRF-14 does not cover is reported, not clamped to the nearest year it does. Both
    // truncations are driven: rotations_at is a template, clang counts its instantiations
    // separately, and the differential tests below use Igrf<10> while everything else uses the
    // default Igrf<13>.
    const ib::Igrf<10> m10 = *ib::Igrf<10>::at(2015.5);
    for (int year : {1850, 2099}) {
        EXPECT_EQ(api::rotations_at(year, 1, 0.0, m).status, Status::DomainError) << year;
        EXPECT_EQ(api::rotations_at(year, 1, 0.0, m10).status, Status::DomainError) << year;
    }
    EXPECT_EQ(api::rotations_at(kYear, kDoy, kUt, m10).status, Status::Ok);
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

    // A non-finite input is reported rather than propagated as a NaN hour — each of the two SM
    // components that enter the clock angle, since `||` short-circuits and one test would leave the
    // other check unexercised.
    const double nan = std::numeric_limits<double>::quiet_NaN();
    for (const fx::vec3d& v : {fx::vec3d{nan, 1.0, 1.0}, fx::vec3d{1.0, nan, 1.0},
                               fx::vec3d{1.0, 1.0, nan},
                               fx::vec3d{std::numeric_limits<double>::infinity(), 1.0, 1.0}}) {
        const Result<double> bad = api::get_mlt(Position<Frame::GEO>{v}, id);
        EXPECT_EQ(bad.status, Status::DomainError);
        EXPECT_EQ(bad.value, 0.0);
    }

    // The other side of midnight. A negative-zero ordinate makes atan2 return -pi rather than +pi,
    // and a negative zero DOES survive the rotation (the dot product sums three of them), so this
    // input is reachable and get_mlt has no lower fold to catch it. What makes that safe is an
    // exact arithmetic fact, asserted here rather than assumed: (-pi) * (12/pi) is exactly -12 in
    // IEEE double, so the -pi end of the range evaluates to exactly +0.0 and never undershoots.
    const Position<Frame::SM> minus_zero =
        ib::transform<Frame::SM>(Position<Frame::GEO>{fx::vec3d{-7.0, -0.0, -0.0}}, id);
    EXPECT_TRUE(std::signbit(minus_zero.v[1]));
    EXPECT_EQ(std::atan2(minus_zero.v[1], minus_zero.v[0]), -std::numbers::pi);
    EXPECT_EQ(12.0 + (-std::numbers::pi * api::hours_per_radian), 0.0);
    const Result<double> other_side =
        api::get_mlt(Position<Frame::GEO>{fx::vec3d{-7.0, -0.0, -0.0}}, id);
    EXPECT_EQ(other_side.status, Status::Ok);
    EXPECT_EQ(other_side.value, 0.0);
    EXPECT_FALSE(std::signbit(other_side.value));
}

TEST(IrbemApi, MltMatchesTheReference) {
    // At MATCHED options — the reference's mid-year IGRF sampling, see rotations_matched — what is
    // left is the two libraries' different solar-ephemeris and sidereal-time series, and nothing
    // else. Measured maximum over the sweep: 9.6e-5 h, i.e. 0.34 s of local time, or 1.4e-3 deg of
    // SM azimuth. The cap is set just above that and the maximum is PRINTED, so a regression that
    // doubles it is visible rather than merely under cap.
    double worst = 0.0;
    double worst_unmatched = 0.0;
    const auto circular_gap = [](double a, double b) {
        const double d = std::abs(a - b);
        return d > 12.0 ? 24.0 - d : d;  // hour angles on a 24 h circle: 23.99 and 0.01 are 0.02 apart
    };
    for (const MltGolden& g : kMltGoldens) {
        const Result<double> mine = api::get_mlt(Position<Frame::GEO>{g.geo},
                                                 rotations_matched(g.year, g.doy, g.ut));
        ASSERT_EQ(mine.status, Status::Ok);
        const double d = circular_gap(mine.value, g.mlt);
        worst = std::max(worst, d);
        EXPECT_LT(d, 2e-4) << g.year << "/" << g.doy << " @" << g.ut;

        // The same comparison at the EXACT epoch, which is what a naive differential test would do.
        const Result<double> unmatched =
            api::get_mlt(Position<Frame::GEO>{g.geo}, rotations_for(g.year, g.doy, g.ut));
        worst_unmatched = std::max(worst_unmatched, circular_gap(unmatched.value, g.mlt));
    }
    std::printf("[ MLT vs oracle ] max |dMLT| = %.4g h (%.4g s) matched, %.4g h unmatched, %zu points\n",
                worst, worst * 3600.0, worst_unmatched, std::size(kMltGoldens));
    EXPECT_GT(worst, 0.0);  // a zero here would mean the goldens were not actually compared
    // The convention difference is real and an order larger than the arithmetic difference: this is
    // the assertion that keeps rotations_matched from looking like an arbitrary fudge.
    EXPECT_GT(worst_unmatched, 10.0 * worst);
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
            if (a == b) { EXPECT_EQ(there.value, in.value); }
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
    // Nine output frames x eight epochs from 1965 to 2029 for one seed, plus three more seeds at
    // three epochs, plus a few non-GEO inputs, all at MATCHED options (see rotations_matched). The
    // tolerance is PER FRAME rather than one blanket number, because each leg has a different and
    // nameable cause of disagreement — a blanket 1e-3 would let a genuine regression in the dipole
    // algebra hide behind the solar ephemeris. Every figure below is the worst over THIS golden
    // set; the number in brackets is the worst over a wider 296-point x 8-epoch sweep (radii
    // 1.02-10 Re, latitudes +/-89 deg, longitudes +/-175 deg) run against the same oracle, quoted
    // so the caps are not mistaken for the routine's whole envelope:
    //
    //   GEO        exact for a GEO input; 3.6e-11 for a GDZ/SPH/RLL input, which is that leg's.
    //   SPH        6.7e-16 [2.4e-16] — bit level, but only once the longitude conventions match.
    //              The reference folds SPH into [0, 360) and leaves GDZ and RLL in (-180, 180];
    //              before api.hpp did the same the raw gap was a full 360 deg for every point with
    //              a negative GEO y. See CoordTransUsesTheReferenceLongitudeConventions.
    //   MAG        5.6e-16 [1.6e-15] — bit level. Hapgood T5 is a pure dipole rotation, so with the
    //              coefficients matched there is nothing left for the two to disagree about. This
    //              is the strongest single piece of evidence that the rotation chain is right.
    //   GDZ, RLL   3.3e-8 and 1.5e-8 [4.2e-8 for RLL]. The reference solves geodetic latitude
    //              iteratively; coords_geodetic.hpp uses the closed form. The gap GROWS towards the
    //              surface — 2.1e-6 km of altitude at 1.02 Re against 1.5e-4 km at 1.5 Re — so the
    //              3e-10 a single geosynchronous seed shows is not the number to quote, and a cap
    //              set from that seed alone is a cap the next low-altitude golden breaks.
    //   GEI        7.2e-5 [1.9e-4]. Sidereal time: IAU 1982 here against the reference's own series.
    //   GSM, SM    5.8e-5 [2.5e-4]. The solar ephemeris, partly cancelled because the dipole axis
    //              is carried through the same T2 that carries the point.
    //   GSE        3.1e-4 [6.3e-4]. The solar ephemeris alone, uncancelled — the largest residual,
    //              and it is a model difference rather than an implementation one.
    constexpr std::array<double, 9> kTolerance{1e-7, 1e-15, 2e-4, 5e-4, 2e-4, 2e-4, 1e-14, 1e-14, 1e-7};

    std::array<double, 9> worst{};
    std::array<double, 9> applied{};  // the cap actually applied, which a non-GEO input can loosen
    for (const CoordGolden& g : kCoordGoldens) {
        const Result<fx::vec3d> mine = api::coord_trans(g.sysaxes_in, g.sysaxes_out, g.in,
                                                        rotations_matched(g.year, g.doy, g.ut));
        ASSERT_EQ(mine.status, Status::Ok);
        const double d = relative_gap(mine.value, g.out);
        const std::size_t out_slot = static_cast<std::size_t>(g.sysaxes_out);
        worst[out_slot] = std::max(worst[out_slot], d);
        // A non-GEO input carries its own leg's error in as well, so the pair is judged by the
        // looser of the two legs.
        const double cap = std::max(kTolerance[static_cast<std::size_t>(g.sysaxes_in)],
                                    kTolerance[out_slot]);
        applied[out_slot] = std::max(applied[out_slot], cap);
        EXPECT_LE(d, cap) << g.year << "/" << g.doy << " " << g.sysaxes_in << "->" << g.sysaxes_out;
    }
    static constexpr std::array<const char*, 9> kNames{"GDZ", "GEO", "GSM", "GSE", "SM",
                                                       "GEI", "MAG", "SPH", "RLL"};
    for (std::size_t i = 0; i < 9; ++i) {
        // `applied` differs from kTolerance for ->GEO only, where the rows with a GDZ/SPH/RLL input
        // carry that leg's iterative-latitude gap in with them.
        std::printf("[ coord_trans vs oracle ] ->%-3s worst = %.4g relative (cap %.4g)\n", kNames[i],
                    worst[i], applied[i]);
    }
    // MAG agreeing to a fraction of an ulp is the load-bearing observation; assert it rather than
    // only printing it, so a change that quietly perturbs the dipole rotation fails here.
    EXPECT_LT(worst[6], 1e-14);
    // SPH likewise: it is bit level only because the longitude fold matches, and if that fold is
    // ever removed this assertion fails at 1e0 rather than drifting quietly.
    EXPECT_LT(worst[7], 1e-14);
    EXPECT_GT(worst[3], 0.0);  // and the GSE leg really was compared
}

TEST(IrbemApi, CoordTransUsesTheReferenceLongitudeConventions) {
    // The reference's OWN two angular conventions disagree, and this test is the only thing in the
    // suite that can see it: every value below was measured on coord_trans_vec1_ at 2015-180
    // 12:00 UT, and the SPH row differs from the RLL and GDZ rows by exactly 360 degrees.
    //
    // Before api.hpp folded SPH to match, this leg was silently 360 deg out for every point with a
    // negative GEO y — half the sky — and the whole golden set could not see it, because it drew
    // its GEO seed from one point whose longitude is +165. Folding SPH to [0, 360) inside
    // detail::sph_in_irbem_longitude is what makes ->SPH agree at 6.7e-16 in
    // CoordTransMatchesTheReference rather than at 1e0.
    const Rotations r = rotations_matched(kYear, kDoy, kUt);
    const fx::vec3d negative_y{1.5, -1.5, 4.0};

    const Result<fx::vec3d> sph = api::coord_trans(1, 7, negative_y, r);
    ASSERT_EQ(sph.status, Status::Ok);
    EXPECT_DOUBLE_EQ(sph.value[0], 4.5276925690687087);
    EXPECT_DOUBLE_EQ(sph.value[1], 62.06164727039765);
    EXPECT_DOUBLE_EQ(sph.value[2], 315.0);  // the reference's SPH answer, NOT -45

    for (int angular : {0, 8}) {
        const Result<fx::vec3d> a = api::coord_trans(1, angular, negative_y, r);
        ASSERT_EQ(a.status, Status::Ok) << angular;
        // GDZ and RLL keep (-180, 180], which is also what the reference returns for codes 0 and 8.
        EXPECT_DOUBLE_EQ(a.value[2], -45.0) << angular;
    }

    // The typed leg is NOT folded: it is this module's own surface and keeps this module's own
    // convention. Naming the difference is the point — the two answers are one meridian.
    const Position<Frame::SPH> typed = api::car2sph(Position<Frame::GEO>{negative_y});
    EXPECT_DOUBLE_EQ(typed.longitude(), -45.0);
    EXPECT_DOUBLE_EQ(typed.longitude() + 360.0, sph.value[2]);

    // The fold is `< 0`, not `<= 0`. Measured: the reference returns -0 for GEO (1, -0, 0) and
    // +180 for (-1, +-0, 0), and both are reproduced by that one compare.
    const Result<fx::vec3d> minus_zero = api::coord_trans(1, 7, fx::vec3d{1.0, -0.0, 0.0}, r);
    ASSERT_EQ(minus_zero.status, Status::Ok);
    EXPECT_EQ(minus_zero.value[2], 0.0);
    EXPECT_TRUE(std::signbit(minus_zero.value[2]));
    for (double y : {0.0, -0.0}) {
        const Result<fx::vec3d> anti = api::coord_trans(1, 7, fx::vec3d{-1.0, y, 0.0}, r);
        ASSERT_EQ(anti.status, Status::Ok);
        EXPECT_DOUBLE_EQ(anti.value[2], 180.0);
    }

    // ...and folding does not break the way back: sph_to_car is periodic, so 315 and -45 land on
    // the same GEO point to the last few bits.
    const Result<fx::vec3d> back = api::coord_trans(7, 1, sph.value, r);
    ASSERT_EQ(back.status, Status::Ok);
    EXPECT_LT(relative_gap(back.value, negative_y), 1e-15);
}

TEST(IrbemApi, CoordTransIsRightWhereTheReferenceIsWrongOnThePolarAxis) {
    // A divergence that is NOT a residual: on the geographic polar axis the reference's geodetic
    // altitude is simply wrong, and this module's is exact. At latitude 90 the geodetic height is
    // |z| - b, with b the WGS84 semi-minor axis, because the ellipsoid normal at the pole IS the
    // axis. Measured on coord_trans_vec1_(1 -> 0) at 2015-180 12:00 UT:
    //
    //     GEO (0, 0, 1.0) Re   reference 0 km            truth 14.447685755... km
    //     GEO (0, 0, 1.5) Re   reference 3178.376157 km  truth 3200.047685755... km
    //     GEO (0, 0, 0)        reference latitude 180    truth latitude 0 (degenerate)
    //
    // 14.4 km and 21.7 km respectively. RLL and SPH are unaffected: both carry a radius rather than
    // an altitude, and both reference legs agree with this module at the pole.
    const Rotations r = rotations_matched(kYear, kDoy, kUt);
    const double b_km = ib::detail::wgs84_semi_minor_km;
    const double re_km = ib::detail::earth_radius_km;

    for (double z : {1.0, 1.5, 3.0, -1.5}) {
        const Result<fx::vec3d> gdz = api::coord_trans(1, 0, fx::vec3d{0.0, 0.0, z}, r);
        ASSERT_EQ(gdz.status, Status::Ok);
        // 1e-9 km is a nanometre; this is `==` in all but name, and it is written as a tolerance
        // only because the km conversion is a multiply the compiler may contract.
        EXPECT_NEAR(gdz.value[0], (std::abs(z) * re_km) - b_km, 1e-9) << z;
        EXPECT_DOUBLE_EQ(gdz.value[1], z > 0.0 ? 90.0 : -90.0) << z;
    }
    // The reference's two answers, frozen, so that "we differ here" stays a measurement.
    EXPECT_GT(std::abs(api::coord_trans(1, 0, fx::vec3d{0.0, 0.0, 1.0}, r).value[0] - 0.0), 14.0);
    EXPECT_GT(std::abs(api::coord_trans(1, 0, fx::vec3d{0.0, 0.0, 1.5}, r).value[0] - 3178.376157),
              21.0);

    // RLL and SPH on the axis: radius and latitude both reproduce the reference exactly.
    for (int code : {7, 8}) {
        const Result<fx::vec3d> a = api::coord_trans(1, code, fx::vec3d{0.0, 0.0, 1.5}, r);
        ASSERT_EQ(a.status, Status::Ok) << code;
        EXPECT_DOUBLE_EQ(a.value[0], 1.5) << code;
        EXPECT_DOUBLE_EQ(a.value[1], 90.0) << code;
    }
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

    // Per-point epochs take the OTHER branch and must agree with the scalar lane just as exactly.
    // Every instantiation is driven down both branches: the hoisted one and the per-point one are
    // different code, and clang counts them per instantiation, so exercising each specialization on
    // only one branch would leave the other unmeasured in a suite that reported 100%.
    std::vector<Rotations> per_point;
    for (std::size_t i = 0; i < geo.size(); ++i) {
        per_point.push_back(rotations_for(1980 + static_cast<int>(i), 100, 3600.0));
    }
    std::vector<Position<Frame::GSM>> gsm_varied(geo.size());
    const Result<bool> a2 = api::transform_vec<Frame::GSM>(
        std::span<const Position<Frame::GEO>>{geo}, std::span<Position<Frame::GSM>>{gsm_varied},
        per_point);
    ASSERT_EQ(a2.status, Status::Ok);
    for (std::size_t i = 0; i < geo.size(); ++i) {
        EXPECT_EQ(gsm_varied[i], ib::transform<Frame::GSM>(geo[i], per_point[i]));
    }
    EXPECT_NE(gsm_varied.front(), gsm.front());  // the epochs really do differ

    std::vector<FieldVector<Frame::MAG>> mag_varied(field.size());
    const Result<bool> b2 = api::transform_vec<Frame::MAG>(
        std::span<const FieldVector<Frame::GEO>>{field},
        std::span<FieldVector<Frame::MAG>>{mag_varied}, per_point);
    ASSERT_EQ(b2.status, Status::Ok);
    for (std::size_t i = 0; i < field.size(); ++i) {
        EXPECT_EQ(mag_varied[i], ib::transform<Frame::MAG>(field[i], per_point[i]));
    }

    std::vector<Position<Frame::GEO>> back(geo.size());
    const Result<bool> c = api::transform_vec<Frame::GEO>(
        std::span<const Position<Frame::GSM>>{gsm}, std::span<Position<Frame::GEO>>{back},
        per_point);
    ASSERT_EQ(c.status, Status::Ok);
    for (std::size_t i = 0; i < gsm.size(); ++i) {
        EXPECT_EQ(back[i], ib::transform<Frame::GEO>(gsm[i], per_point[i]));
    }
    const Result<bool> c2 = api::transform_vec<Frame::GEO>(
        std::span<const Position<Frame::GSM>>{gsm}, std::span<Position<Frame::GEO>>{back}, shared);
    ASSERT_EQ(c2.status, Status::Ok);
    for (std::size_t i = 0; i < gsm.size(); ++i) {
        // GEO -> GSM -> GEO at one epoch is the identity to a rounding, which is the round trip the
        // hoisted path must not disturb.
        EXPECT_EQ(back[i], ib::transform<Frame::GEO>(gsm[i], r));
        EXPECT_LT(relative_gap(back[i].v, geo[i].v), 4e-15);
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

    // Every transform_vec instantiation the suite creates gets both rejections: clang counts
    // template instantiations separately, so checking one specialization would leave the others'
    // guards unmeasured behind a report that said 100%.
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

    std::vector<Position<Frame::GEO>> geo(4);
    std::vector<Position<Frame::GEO>> short_geo(3);
    EXPECT_EQ((api::transform_vec<Frame::GEO>(std::span<const Position<Frame::GSM>>{gsm},
                                              std::span<Position<Frame::GEO>>{short_geo}, shared)
                   .status),
              Status::DomainError);
    EXPECT_EQ((api::transform_vec<Frame::GEO>(std::span<const Position<Frame::GSM>>{gsm},
                                              std::span<Position<Frame::GEO>>{geo}, two)
                   .status),
              Status::DomainError);

    std::vector<FieldVector<Frame::GEO>> b_geo(4);
    std::vector<FieldVector<Frame::MAG>> b_mag(4);
    std::vector<FieldVector<Frame::MAG>> short_mag(3);
    EXPECT_EQ((api::transform_vec<Frame::MAG>(std::span<const FieldVector<Frame::GEO>>{b_geo},
                                              std::span<FieldVector<Frame::MAG>>{short_mag}, shared)
                   .status),
              Status::DomainError);
    EXPECT_EQ((api::transform_vec<Frame::MAG>(std::span<const FieldVector<Frame::GEO>>{b_geo},
                                              std::span<FieldVector<Frame::MAG>>{b_mag}, two)
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

// ---- the MAKE_LSTAR family -------------------------------------------------------------------

TEST(IrbemApi, LstarPhiInvertsItself) {
    const ib::Igrf<> m = *ib::Igrf<>::at(2015.5);
    const double k0 = ib::dipole_moment(m);

    // Phi = 2*pi*k0/L*, so the two directions are the same expression and must compose to identity.
    for (double lstar : {1.5, 4.0, 6.6, 12.0}) {
        const Result<double> phi = api::lstar_phi(1, lstar, k0);
        ASSERT_EQ(phi.status, Status::Ok);
        const Result<double> back = api::lstar_phi(2, phi.value, k0);
        ASSERT_EQ(back.status, Status::Ok);
        EXPECT_NEAR(back.value, lstar, 1e-12 * lstar);
        EXPECT_GT(phi.value, 0.0);
    }
    // Larger L* means a SMALLER enclosed flux — the relation is a reciprocal, and getting its sense
    // backwards is the one way to be wrong here that still produces plausible numbers.
    EXPECT_LT(api::lstar_phi(1, 6.6, k0).value, api::lstar_phi(1, 4.0, k0).value);

    // The direction code is checked, not defaulted, exactly as `sysaxes` is.
    for (int bad : {-1, 0, 3, 99}) EXPECT_EQ(api::lstar_phi(bad, 4.0, k0).status, Status::DomainError);
    // Zero has no image under a reciprocal, and neither does a negative moment.
    EXPECT_EQ(api::lstar_phi(1, 0.0, k0).status, Status::DomainError);
    EXPECT_EQ(api::lstar_phi(1, -4.0, k0).status, Status::DomainError);
    EXPECT_EQ(api::lstar_phi(1, std::numeric_limits<double>::infinity(), k0).status,
              Status::DomainError);
    EXPECT_EQ(api::lstar_phi(1, 4.0, 0.0).status, Status::DomainError);
    // NaN and +inf reach the moment guard by different routes: `!(NaN > 0)` is true and
    // short-circuits, while `+inf > 0` is true and only the isfinite test catches it.
    EXPECT_EQ(api::lstar_phi(1, 4.0, std::numeric_limits<double>::quiet_NaN()).status,
              Status::DomainError);
    EXPECT_EQ(api::lstar_phi(1, 4.0, std::numeric_limits<double>::infinity()).status,
              Status::DomainError);
}

TEST(IrbemApi, LstarPhiMatchesTheReference) {
    // LstarPhiInvertsItself cannot see the CONSTANT: self-inversion, the reciprocal ordering, the
    // sign and every domain guard all survive replacing 2*pi with pi, or with 1, or with 1e9.
    // Verified by perturbation — dropping the 2 from `2.0 * pi` leaves all of that green. So the
    // absolute value has to be pinned against something outside this file, and the reference is
    // right there: lstar_phi1_(ntime=1, whichinv, options={1,0,0,0,0}, iyear, idoy=180, ...) on
    // libirbem-O2.so, frozen at %.17g.
    //
    // This also settles a question the header used to answer wrongly. The reference takes iyear and
    // idoy, and the implied B0 below — Phi * L* / 2pi — is 30939.3888 nT in 1965 falling to
    // 29660.6279 in 2029, i.e. it tracks the epoch's own dipole exactly as dipole_moment() does.
    // It is NOT frozen at a 1960s value, and there is no convention here to reconcile.
    struct PhiGolden {
        int year;    ///< `iyear`.
        double phi;  ///< `Phi` for `Lstar = 4`, nT Re^2, as the reference returned it.
    };
    constexpr PhiGolden kPhiGoldens[] = {
        {1965, 48599.478358059241},
        {1990, 47607.432601650558},
        {2015, 46905.583704727214},
        {2029, 46590.805305975191},
    };
    for (const PhiGolden& g : kPhiGoldens) {
        const std::optional<ib::Igrf<>> m = ib::Igrf<>::at(static_cast<double>(g.year) + 0.5);
        ASSERT_TRUE(m.has_value()) << g.year;
        const double k0 = ib::dipole_moment(*m);
        const Result<double> phi = api::lstar_phi(1, 4.0, k0);
        ASSERT_EQ(phi.status, Status::Ok) << g.year;
        // 1.5e-16 is the worst of the four, i.e. one ulp. A blanket relative tolerance is used
        // rather than `==` only because the reference's own multiply order is not ours.
        EXPECT_LT(std::abs(phi.value - g.phi) / g.phi, 4e-16) << g.year;
        // ...and the reverse direction, which the reference also serves.
        const Result<double> back = api::lstar_phi(2, g.phi, k0);
        ASSERT_EQ(back.status, Status::Ok) << g.year;
        EXPECT_LT(std::abs(back.value - 4.0), 1e-15) << g.year;
    }
}

TEST(IrbemApi, MakeLstarVecDeviceLaneAgreesWithTheHostLane) {
    // The question MakeLstarVecAgreesWithTheScalarLane deliberately does NOT ask, because it pins
    // the host lane to compare with `==`. Here the auto lane runs as production runs it and is held
    // to docs/ERROR_BUDGET.md instead. On a machine with no device the two runs are the same lane
    // and the assertions are trivially satisfied; the printout says which happened, so a green run
    // never implies a device was exercised.
    const ib::Igrf<> m = *ib::Igrf<>::at(2015.5);
    const Rotations r = rotations_matched(kYear, kDoy, kUt);
    std::vector<Position<Frame::GEO>> pts;
    for (int i = 0; i < 16; ++i) {
        const double t = static_cast<double>(i) / 15.0;
        const double l = 3.0 + (2.0 * t);
        const double ph = 2.0 * std::numbers::pi * t * 3.0;
        pts.push_back(Position<Frame::GEO>{fx::vec3d{l * std::cos(ph), l * std::sin(ph), 0.3}});
    }
    std::vector<api::MagneticCoordinates> automatic(pts.size()), host(pts.size());
    std::vector<api::MagneticCoordinates> again(pts.size());
    std::vector<Status> sa(pts.size()), sh(pts.size()), s2(pts.size());

    const Result<bool> ra = api::make_lstar_vec(m, r, pts, automatic, sa);
    const Result<bool> r2 = api::make_lstar_vec(m, r, pts, again, s2);
    ASSERT_EQ(ra.status, Status::Ok);
    // Whichever lane ran, it must be REPRODUCIBLE: the flux reduction is ordered and in fp64 by
    // ERROR_BUDGET.md 6, so two identical calls must return identical bits even on the device.
    for (std::size_t i = 0; i < pts.size(); ++i) EXPECT_EQ(automatic[i], again[i]) << i;
    EXPECT_EQ(ra.value, r2.value);

    double worst_lstar = 0.0;
    double worst_lm = 0.0;
    double worst_b = 0.0;
    {
        const HostLanePin host_only;
        const Result<bool> rh = api::make_lstar_vec(m, r, pts, host, sh);
        ASSERT_EQ(rh.status, Status::Ok);
        // The host lane is the fp64 reference and never reports a device.
        EXPECT_FALSE(rh.value);
    }
    for (std::size_t i = 0; i < pts.size(); ++i) {
        ASSERT_EQ(sa[i], sh[i]) << i;
        if (sa[i] != Status::Ok) continue;
        worst_lstar = std::max(worst_lstar, std::abs(automatic[i].lstar - host[i].lstar));
        worst_lm = std::max(worst_lm, std::abs(automatic[i].lm - host[i].lm));
        worst_b = std::max(worst_b,
                           std::abs(automatic[i].blocal - host[i].blocal) / host[i].blocal);
        // MLT is pure host geometry in both lanes, so it must be identical to the bit.
        EXPECT_EQ(automatic[i].mlt, host[i].mlt) << i;
    }
    std::printf("[ make_lstar_vec lanes ] auto lane served by device: %s ; max |dL*| = %.3g, "
                "|dLm| = %.3g, rel |dB| = %.3g over %zu points\n",
                ra.value ? "YES" : "no", worst_lstar, worst_lm, worst_b, pts.size());
    // ERROR_BUDGET.md 4: 0.01 absolute on L*. Lm's own row is 1e-3, and the device lane does NOT
    // always hold it — measured 4.0e-3 over 512 points spread across L = 3..5 — so the cap here is
    // the L* one for both, and the gap is REPORTED rather than hidden. A caller who needs Lm to
    // 1e-3 pins the host lane; api.hpp's make_lstar_vec brief says so.
    EXPECT_LT(worst_lstar, 0.01);
    EXPECT_LT(worst_lm, 0.01);
    // B is evaluated, not integrated, so it is held to the Blocal row: 1e-6 relative.
    EXPECT_LT(worst_b, 1e-6);
}

TEST(IrbemApi, MakeLstarReportsAnUnsupportedExternalModel) {
    // The whole point of this test: an internal-only L* at geosynchronous during a storm is not a
    // slightly worse T89 L*, it is a different number. Asking for an external field this module
    // cannot supply must be REFUSED, never silently answered with the internal field.
    const ib::Igrf<> m = *ib::Igrf<>::at(2015.5);
    const Rotations r = rotations_matched(kYear, kDoy, kUt);
    const std::array<Position<Frame::GEO>, 2> pts{Position<Frame::GEO>{fx::vec3d{4.0, 0.0, 0.0}},
                                                  Position<Frame::GEO>{fx::vec3d{0.0, 5.0, 0.0}}};
    std::array<api::MagneticCoordinates, 2> out{};
    std::array<Status, 2> st{};

    // Igrf<10> too: the differential test below instantiates the L* lanes at that truncation, and a
    // guard is only measured in the instantiation it is compiled into.
    const ib::Igrf<10> m10 = *ib::Igrf<10>::at(2015.5);
    std::array<api::MagneticCoordinates, 2> out10{};
    std::array<Status, 2> st10{};
    EXPECT_EQ(api::make_lstar_vec(m10, r, pts, out10, st10, ib::ExternalModel::Tsyganenko1989).status,
              Status::ParametersMissing);
    EXPECT_TRUE(ib::is_baddata(out10[0].lstar));
    EXPECT_EQ(api::make_lstar_vec(m10, r, pts, out10, st10, static_cast<ib::ExternalModel>(15)).status,
              Status::DomainError);
    EXPECT_EQ(api::make_lstar(m10, r, pts[0], ib::ExternalModel::Tsyganenko1996).status,
              Status::ParametersMissing);
    EXPECT_EQ(api::make_lstar(m10, r, pts[0], static_cast<ib::ExternalModel>(15)).status,
              Status::DomainError);
    std::array<api::MagneticCoordinates, 1> short10{};
    std::array<Status, 1> short_st10{};
    EXPECT_EQ(api::make_lstar_vec(m10, r, pts, short10, st10).status, Status::DomainError);
    EXPECT_EQ(api::make_lstar_vec(m10, r, pts, out10, short_st10).status, Status::DomainError);
    {
        std::span<const Position<Frame::GEO>> empty;
        std::span<api::MagneticCoordinates> empty_out;
        std::span<Status> empty_st;
        EXPECT_EQ(api::make_lstar_vec(m10, r, empty, empty_out, empty_st).status, Status::Ok);
    }

    for (ib::ExternalModel kext : {ib::ExternalModel::Tsyganenko1989, ib::ExternalModel::Tsyganenko1996,
                                   ib::ExternalModel::Tsyganenko2004Storm,
                                   ib::ExternalModel::OlsonPfitzerQuiet1977}) {
        const Result<bool> v = api::make_lstar_vec(m, r, pts, out, st, kext);
        EXPECT_EQ(v.status, Status::ParametersMissing) << ib::name_of(kext);
        EXPECT_EQ(st[0], Status::ParametersMissing);
        // Refused means NOTHING computed — every output is the sentinel, not a plausible number a
        // caller could mistake for an answer under that model's name.
        EXPECT_TRUE(ib::is_baddata(out[0].lstar));
        EXPECT_TRUE(ib::is_baddata(out[0].lm));
        EXPECT_TRUE(ib::is_baddata(out[1].mlt));

        const Result<api::MagneticCoordinates> one = api::make_lstar(m, r, pts[0], kext);
        EXPECT_EQ(one.status, Status::ParametersMissing);
        EXPECT_TRUE(ib::is_baddata(one.value.lstar));
    }
    // A key outside IRBEM's own 0..14 is a different defect and gets a different status.
    EXPECT_EQ(api::make_lstar(m, r, pts[0], static_cast<ib::ExternalModel>(200)).status,
              Status::DomainError);
    EXPECT_EQ(api::make_lstar_vec(m, r, pts, out, st, static_cast<ib::ExternalModel>(15)).status,
              Status::DomainError);

    // Shell splitting refuses on the same terms.
    const std::array<double, 2> alpha{45.0, 90.0};
    std::array<api::MagneticCoordinates, 4> grid{};
    std::array<Status, 4> grid_st{};
    EXPECT_EQ(api::make_lstar_shell_splitting(m, r, pts, alpha, grid, grid_st,
                                              ib::ExternalModel::Tsyganenko1989)
                  .status,
              Status::ParametersMissing);
    EXPECT_TRUE(ib::is_baddata(grid[3].lstar));

    // ...and a length mismatch is caught before any of that.
    std::array<api::MagneticCoordinates, 1> short_out{};
    std::array<Status, 1> short_st{};
    EXPECT_EQ(api::make_lstar_vec(m, r, pts, short_out, st).status, Status::DomainError);
    EXPECT_EQ(api::make_lstar_vec(m, r, pts, out, short_st).status, Status::DomainError);
    EXPECT_EQ(api::make_lstar_shell_splitting(m, r, pts, alpha, short_out, grid_st).status,
              Status::DomainError);
    EXPECT_EQ(api::make_lstar_shell_splitting(m, r, pts, alpha, grid, short_st).status,
              Status::DomainError);

    // An empty batch is a valid batch.
    std::span<const Position<Frame::GEO>> none;
    std::span<api::MagneticCoordinates> no_out;
    std::span<Status> no_st;
    EXPECT_EQ(api::make_lstar_vec(m, r, none, no_out, no_st).status, Status::Ok);
    EXPECT_EQ(api::make_lstar_shell_splitting(m, r, none, alpha, no_out, no_st).status, Status::Ok);
}

TEST(IrbemApi, MakeLstarVecAgreesWithTheScalarLane) {
    // Both lanes on the host: the `==` below is a claim about the ALGEBRA, not about which device
    // happened to service which batch. See HostLanePin.
    const HostLanePin host_only;
    const ib::Igrf<> m = *ib::Igrf<>::at(2015.5);
    const Rotations r = rotations_matched(kYear, kDoy, kUt);
    std::vector<Position<Frame::GEO>> pts;
    for (int i = 0; i < 6; ++i) {
        const double x = 3.0 + (0.5 * static_cast<double>(i));
        pts.push_back(Position<Frame::GEO>{fx::vec3d{x, 0.7, -0.4}});
    }
    std::vector<api::MagneticCoordinates> out(pts.size());
    std::vector<Status> st(pts.size());

    const Result<bool> batch = api::make_lstar_vec(m, r, pts, out, st);
    ASSERT_EQ(batch.status, Status::Ok);
    for (std::size_t i = 0; i < pts.size(); ++i) {
        // `==` on the whole record: the batch and the scalar lane run the same root-find on the
        // same shell, so a tolerance here would hide exactly the reordering bug it looks for.
        const Result<api::MagneticCoordinates> one = api::make_lstar(m, r, pts[i]);
        EXPECT_EQ(one.status, st[i]);
        EXPECT_EQ(one.value, out[i]) << i;
        // ...and the six fields really are the six IRBEM names, in the right slots.
        EXPECT_GT(out[i].lstar, 0.0);
        EXPECT_GE(out[i].blocal, out[i].bmin);      // the local field is never below the minimum
        // L_m and L* are within a few percent of each other for an internal-only field, but NOT
        // ordered: the reference's own goldens have L_m above L* at four points and below it at a
        // fifth, so an inequality here would be asserting a rule that does not exist.
        EXPECT_LT(std::abs(out[i].lm - out[i].lstar) / out[i].lm, 0.05);
        EXPECT_GE(out[i].xj, 0.0);
        EXPECT_EQ(out[i].mlt, api::get_mlt(pts[i], r).value);
    }
}

TEST(IrbemApi, ShellSplittingVariesLstarWithPitchAngle) {
    // The 90-degree column is compared to make_lstar with `==`, so both lanes must be the host's.
    const HostLanePin host_only;
    const ib::Igrf<> m = *ib::Igrf<>::at(2015.5);
    const Rotations r = rotations_matched(kYear, kDoy, kUt);
    const std::array<Position<Frame::GEO>, 2> pts{Position<Frame::GEO>{fx::vec3d{4.0, 0.5, 0.3}},
                                                  Position<Frame::GEO>{fx::vec3d{0.4, 5.0, -0.6}}};
    const std::array<double, 3> alpha{40.0, 65.0, 90.0};
    std::array<api::MagneticCoordinates, 6> grid{};
    std::array<Status, 6> st{};

    const Result<bool> v = api::make_lstar_shell_splitting(m, r, pts, alpha, grid, st);
    ASSERT_EQ(v.status, Status::Ok);

    for (std::size_t i = 0; i < pts.size(); ++i) {
        for (std::size_t k = 0; k < alpha.size(); ++k) {
            const api::MagneticCoordinates& c = grid[(i * alpha.size()) + k];
            // Point-major indexing, as documented: the 90-degree column must reproduce make_lstar,
            // which is what pins the layout down rather than merely describing it.
            if (alpha[k] == 90.0) { EXPECT_EQ(c, api::make_lstar(m, r, pts[i]).value); }
            EXPECT_GT(c.lstar, 0.0);
            // B_local and B_min are properties of the LINE, not of the particle, so they must not
            // move with pitch angle. B_mirror and I must.
            EXPECT_EQ(c.blocal, grid[i * alpha.size()].blocal);
            EXPECT_EQ(c.bmin, grid[i * alpha.size()].bmin);
            EXPECT_EQ(c.mlt, grid[i * alpha.size()].mlt);
        }
        // Shell splitting itself: a smaller pitch angle mirrors further down the line and drifts on
        // a different shell. If L* came out identical across alpha the routine would be computing
        // the 90-degree shell three times.
        EXPECT_NE(grid[i * alpha.size()].lstar, grid[(i * alpha.size()) + 2].lstar);
        EXPECT_GT(grid[i * alpha.size()].xj, grid[(i * alpha.size()) + 2].xj);
    }
}

TEST(IrbemApi, MakeLstarMatchesTheOracle) {
    // The adapter's job here is to put six numbers in six named slots; driftshell.hpp's own suite
    // owns the convergence question. So this is a SLOT test with a physics tolerance rather than a
    // convergence study: a field swap, a unit slip or a transposed index moves one of these by
    // orders of magnitude, which is what the per-field caps below are sized to catch. They are the
    // ERROR_BUDGET.md figures for matched options — L* 0.01 absolute, L_m 1e-3, XJ 1e-4 relative,
    // B 1e-6 relative — loosened only where measurement says the two libraries genuinely differ.
    double worst_lstar = 0.0;
    double worst_lm = 0.0;
    double worst_blocal = 0.0;
    double worst_bmin = 0.0;
    double worst_xj = 0.0;
    double worst_mlt = 0.0;
    int compared = 0;
    for (const LstarGolden& g : kLstarGoldens) {
        // Igrf<10>, not the default Igrf<13>: IRBEM truncates its internal field at degree 10
        // (igrf.hpp says so, and measurement confirms it). Matched options means matched
        // TRUNCATION too. Measured on `Blocal`, which is a pure field evaluation with no tracing
        // in it, the choice is worth two orders and then some — degree 10 at mid-year reproduces
        // the reference to 5e-16 relative, degree 13 at mid-year to 2e-9, and degree 13 at the
        // EXACT date only to 3e-4. Blocal is therefore the sharpest instrument in this test, and
        // its cap below is set accordingly.
        const ib::Igrf<10> m = *ib::Igrf<10>::at(static_cast<double>(g.year) + 0.5);
        const Result<Rotations> rr = api::rotations_at(g.year, g.doy, g.ut, m);
        ASSERT_EQ(rr.status, Status::Ok);
        const Rotations& r = rr.value;
        const Result<api::MagneticCoordinates> mine =
            api::make_lstar(m, r, Position<Frame::GEO>{g.geo});

        // The reference refuses one of these shells — measured, `make_lstar1_` returns -1e30 (NOT
        // the -1e31 its own documentation names) at the L=10 point. There is nothing to compare
        // against there, and this module converging where the reference gave up is not a defect, so
        // the row is reported and skipped rather than turned into an assertion either way.
        if (g.out.lstar < -1e29) {
            std::printf("[ make_lstar vs oracle ] reference refused L=%.2f at %d/%d; ours: %s "
                        "L*=%.4f\n",
                        fx::norm(g.geo), g.year, g.doy, ib::describe(mine.status).data(),
                        mine.value.lstar);
            continue;
        }
        ASSERT_EQ(mine.status, Status::Ok) << g.year << "/" << g.doy;
        ++compared;
        worst_lstar = std::max(worst_lstar, std::abs(mine.value.lstar - g.out.lstar));
        worst_lm = std::max(worst_lm, std::abs(mine.value.lm - g.out.lm) / g.out.lm);
        worst_xj = std::max(worst_xj, std::abs(mine.value.xj - g.out.xj) / std::max(g.out.xj, 1e-3));
        worst_blocal =
            std::max(worst_blocal, std::abs(mine.value.blocal - g.out.blocal) / g.out.blocal);
        worst_bmin = std::max(worst_bmin, std::abs(mine.value.bmin - g.out.bmin) / g.out.bmin);
        double dmlt = std::abs(mine.value.mlt - g.out.mlt);
        if (dmlt > 12.0) dmlt = 24.0 - dmlt;
        worst_mlt = std::max(worst_mlt, dmlt);

        // Blocal is the one output with no algorithm between the model and the answer, so it is
        // held to the model: a slot swap or a unit slip cannot survive 1e-12.
        EXPECT_LT(std::abs(mine.value.blocal - g.out.blocal) / g.out.blocal, 1e-12) << "Blocal";
        // Bmin, XJ, Lm and L* all come out of the tracer and the root-find, whose step sizes the
        // two libraries do not share. These caps are sized to catch a wrong SLOT (orders of
        // magnitude) rather than to certify convergence — driftshell.hpp's own differential suite
        // owns that question, and ERROR_BUDGET.md's 0.01 absolute on L* is its bar, not this
        // file's. Measured here: L* to 0.036 absolute, Lm to 4.6e-3 relative, Bmin to 1.9e-4,
        // XJ to 0.079 relative.
        EXPECT_LT(std::abs(mine.value.lstar - g.out.lstar), 0.05) << "Lstar";
        EXPECT_LT(std::abs(mine.value.lm - g.out.lm) / g.out.lm, 0.02) << "Lm";
        EXPECT_LT(std::abs(mine.value.bmin - g.out.bmin) / g.out.bmin, 5e-3) << "Bmin";
        EXPECT_LT(std::abs(mine.value.xj - g.out.xj) / std::max(g.out.xj, 1e-3), 0.2) << "XJ";
        EXPECT_LT(dmlt, 2e-4) << "MLT";
    }
    std::printf("[ make_lstar vs oracle ] %d points: |dL*| <= %.3g abs, |dLm|/Lm <= %.3g, "
                "|dBlocal|/B <= %.3g, |dBmin|/B <= %.3g, |dXJ|/XJ <= %.3g, |dMLT| <= %.3g h\n",
                compared, worst_lstar, worst_lm, worst_blocal, worst_bmin, worst_xj, worst_mlt);
    EXPECT_GE(compared, 6);
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
