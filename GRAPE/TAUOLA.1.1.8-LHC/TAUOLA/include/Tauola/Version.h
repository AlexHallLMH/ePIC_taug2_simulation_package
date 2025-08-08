/**
 * This file contains tauola version definitions
 *
 * @date January 21 2020
 */
#ifndef TAUOLA_VERSION_H
#define TAUOLA_VERSION_H

#include <string>

/** @brief TAUOLA version string */
#define TAUOLA_VERSION "1.1.8"

/** @brief TAUOLA version code */
#define TAUOLA_VERSION_CODE 10108

namespace Tauolapp {

/* @brief Get TAUOLA library version string */
inline std::string version() {
    return TAUOLA_VERSION;
}

} // namespace Tauolapp

#endif
