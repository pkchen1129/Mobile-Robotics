#include <gtsam/sam/BearingRangeFactor.h>
#include <gtsam/slam/BetweenFactor.h>
#include <gtsam/slam/dataset.h>
#include <gtsam/geometry/Point3.h>
#include <gtsam/geometry/Pose2.h>
#include <gtsam/geometry/Rot3.h>
#include <gtsam/inference/FactorGraph.h>
#include <gtsam/inference/Symbol.h>
#include <gtsam/nonlinear/NonlinearFactor.h>
#include <gtsam/nonlinear/Values.h>
#include <gtsam/nonlinear/Values-inl.h>
#include <gtsam/linear/Sampler.h>
#include <gtsam/base/GenericValue.h>
#include <gtsam/base/Lie.h>
#include <gtsam/base/Matrix.h>
#include <gtsam/base/types.h>
#include <gtsam/base/Value.h>
#include <gtsam/base/Vector.h>

#include <boost/assign/list_inserter.hpp>
#include <boost/filesystem/operations.hpp>
#include <boost/filesystem/path.hpp>

#include <cmath>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <vector>
#include <string>
using namespace std;
using namespace gtsam;
namespace fs = boost::filesystem;

typedef std::pair<NonlinearFactorGraph::shared_ptr, Values::shared_ptr> GraphAndValues;
#define LINESIZE 81920
#define DATASET_DIR "@CMAKE_CURRENT_SOURCE_DIR@/dataset"
static SharedNoiseModel readNoiseModel(ifstream& is, bool smart,
    NoiseFormat noiseFormat, KernelFunctionType kernelFunctionType) {
  double v1, v2, v3, v4, v5, v6;
  is >> v1 >> v2 >> v3 >> v4 >> v5 >> v6;

  if (noiseFormat == NoiseFormatAUTO) {
    // Try to guess covariance matrix layout
    if (v1 != 0.0 && v2 == 0.0 && v3 != 0.0 && v4 != 0.0 && v5 == 0.0
        && v6 == 0.0) {
      // NoiseFormatGRAPH
      noiseFormat = NoiseFormatGRAPH;
    } else if (v1 != 0.0 && v2 == 0.0 && v3 == 0.0 && v4 != 0.0 && v5 == 0.0
        && v6 != 0.0) {
      // NoiseFormatCOV
      noiseFormat = NoiseFormatCOV;
    } else {
      throw std::invalid_argument(
          "load2D: unrecognized covariance matrix format in dataset file. Please specify the noise format.");
    }
  }

  // Read matrix and check that diagonal entries are non-zero
  Matrix M(3, 3);
  switch (noiseFormat) {
  case NoiseFormatG2O:
  case NoiseFormatCOV:
    // i.e., [ v1 v2 v3; v2' v4 v5; v3' v5' v6 ]
    if (v1 == 0.0 || v4 == 0.0 || v6 == 0.0)
      throw runtime_error(
          "load2D::readNoiseModel looks like this is not G2O matrix order");
    M << v1, v2, v3, v2, v4, v5, v3, v5, v6;
    break;
  case NoiseFormatTORO:
  case NoiseFormatGRAPH:
    // http://www.openslam.org/toro.html
    // inf_ff inf_fs inf_ss inf_rr inf_fr inf_sr
    // i.e., [ v1 v2 v5; v2' v3 v6; v5' v6' v4 ]
    if (v1 == 0.0 || v3 == 0.0 || v4 == 0.0)
      throw invalid_argument(
          "load2D::readNoiseModel looks like this is not TORO matrix order");
    M << v1, v2, v5, v2, v3, v6, v5, v6, v4;
    break;
  default:
    throw runtime_error("load2D: invalid noise format");
  }

  // Now, create a Gaussian noise model
  // The smart flag will try to detect a simpler model, e.g., unit
  SharedNoiseModel model;
  switch (noiseFormat) {
  case NoiseFormatG2O:
  case NoiseFormatTORO:
    // In both cases, what is stored in file is the information matrix
    model = noiseModel::Gaussian::Information(M, smart);
    break;
  case NoiseFormatGRAPH:
  case NoiseFormatCOV:
    // These cases expect covariance matrix
    model = noiseModel::Gaussian::Covariance(M, smart);
    break;
  default:
    throw invalid_argument("load2D: invalid noise format");
  }

  switch (kernelFunctionType) {
  case KernelFunctionTypeNONE:
    return model;
    break;
  case KernelFunctionTypeHUBER:
    return noiseModel::Robust::Create(
        noiseModel::mEstimator::Huber::Create(1.345), model);
    break;
  case KernelFunctionTypeTUKEY:
    return noiseModel::Robust::Create(
        noiseModel::mEstimator::Tukey::Create(4.6851), model);
    break;
  default:
    throw invalid_argument("load2D: invalid kernel function type");
  }
}