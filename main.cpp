#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <vector>

using namespace std;

struct LonLat {
    double lon;
    double lat;
};

ostream& operator<<(ostream& os, const LonLat& p)
{
    os << fixed << setprecision(12) << "{lon[" << p.lon << "] lat[" << p.lat << "]}";
    return os;
}

LonLat ToRad(const LonLat& p)
{
    static constexpr double kDeg2Rad = M_PI / 180.0;
    return { p.lon * kDeg2Rad, p.lat * kDeg2Rad };
}

struct Vec {
    double x;
    double y;
};

ostream& operator<<(ostream& os, const Vec& p)
{
    os << fixed << setprecision(12) << "{x[" << p.x << "] y[" << p.y << "]}";
    return os;
}

static constexpr double kEarthRadiusMean = 6371000.0;
static constexpr double kEarthRadiusMax = 6378137.0;
static constexpr double kEarthE = 0.08181919108425; // 子午圈偏心率
static constexpr double kEarthE2 = kEarthE * kEarthE;

Vec MapProject1(const LonLat& p0, const LonLat& p1)
{
    double sin_lat0 = sin(p0.lat);
    double cos_lat0 = cos(p0.lat);
    double sin_lat1 = sin(p1.lat);
    double cos_lat1 = cos(p1.lat);
    double cos_d_lon = cos(p1.lon - p0.lon);

    double arg = sin_lat0 * sin_lat1 + cos_lat0 * cos_lat1 * cos_d_lon;
    arg = clamp(arg, -1.0, 1.0);
    double c = acos(arg);
    double k = (fabs(c) < numeric_limits<double>::epsilon()) ? 1.0 : (c / sin(c));

    double x = k * cos_lat1 * sin(p1.lon - p0.lon) * kEarthRadiusMean;
    double y = k * (cos_lat0 * sin_lat1 - sin_lat0 * cos_lat1 * cos_d_lon) * kEarthRadiusMean;
    return { x, y };
}

Vec MapProject2(const LonLat& p0, const LonLat& p1)
{
    // this is WRONG!
    double sin_lat1 = sin(p1.lat);
    double sin_lat1_2 = sin_lat1 * sin_lat1;
    double tmp = 1.0 / sqrt(1.0 - kEarthE2 * sin_lat1_2);
    double rn = kEarthRadiusMax * tmp; // 当前点地心距离
    double rm = kEarthRadiusMax * (1 - kEarthE2) * tmp * tmp * tmp; // 子午圈曲率半径
    cout << "\33[31m" << fixed << setprecision(12) << tmp << " " << rn << " " << rm << "\33[0m\n";

    double x = (p1.lon - p0.lon) * rn * cos(p1.lat);
    double y = (p1.lat - p0.lat) * rm;
    return { x, y };
}

Vec MapProject3(const LonLat& p0, const LonLat& p1)
{
    // floatify MapProject2
    float sin_lat1 = sin(float(p1.lat));
    float sin_lat1_2 = sin_lat1 * sin_lat1;
    float tmp = 1.0f / sqrt(1.0 - float(kEarthE2) * sin_lat1_2);
    float rn = float(kEarthRadiusMax) * tmp;
    float rm = float(kEarthRadiusMax) * (1.0f - float(kEarthE2)) * tmp * tmp * tmp;
    double x = float(p1.lon - p0.lon) * cos(float(p1.lat)) * rn;
    double y = float(p1.lat - p0.lat) * rm;
    return { x, y };
}

Vec MapProject4(const LonLat& p0, const LonLat& p1)
{
    // simplification of MapProject2, but correct
    double x = float(p1.lon - p0.lon) * cos(float(p1.lat)) * float(kEarthRadiusMean);
    double y = float(p1.lat - p0.lat) * float(kEarthRadiusMean);
    return { x, y };
}

Vec Diff(const Vec& ref, const Vec& v)
{
    return { (v.x - ref.x) / ref.x, (v.y - ref.y) / ref.y };
}

int main()
{
    vector<LonLat> pts = {
        { 116.280376123, 40.051177456 },
        { 116.280386123, 40.051177456 },
        { 116.280489123, 40.051643456 },
        { 116.233323123, 40.149726456 },
        { 113.331021123, 23.118705456 },
        { 113.331094123, 23.118722456 },
        { 113.355309123, 23.152298456 }
    };
    for (size_t i = 0; i + 1 < pts.size(); ++i) {
        const auto& p1 = pts[i];
        const auto& p2 = pts[i + 1];
        auto p1_rad = ToRad(p1);
        auto p2_rad = ToRad(p2);
        auto r1 = MapProject1(p1_rad, p2_rad);
        auto r2 = MapProject2(p1_rad, p2_rad);
        auto r3 = MapProject3(p1_rad, p2_rad);
        auto r4 = MapProject4(p1_rad, p2_rad);
        auto& ref = r1;
        cout << "p1: " << p1 << "\n"
             << "p2: " << p2 << "\n"
             << "r1: " << r1 << " " << Diff(ref, r1) << "\n"
             << "r2: " << r2 << " " << Diff(ref, r2) << "\n"
             << "r3: " << r3 << " " << Diff(ref, r3) << "\n"
             << "r4: " << r4 << " " << Diff(ref, r4) << "\n\n";
    }
}