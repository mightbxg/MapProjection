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
    os << fixed << setprecision(16) << "{lon[" << p.lon << "] lat[" << p.lat << "]}";
    return os;
}

LonLat ToRad(const LonLat& p)
{
    static constexpr double kDeg2Rad = M_PI / 180.0;
    return { p.lon * kDeg2Rad, p.lat * kDeg2Rad };
}

LonLat ToDeg(const LonLat& p)
{
    static constexpr double kRad2Deg = 180.0 / M_PI;
    return { p.lon * kRad2Deg, p.lat * kRad2Deg };
}

struct Vec {
    double x;
    double y;
};

ostream& operator<<(ostream& os, const Vec& p)
{
    os << fixed << setprecision(16) << "{x[" << p.x << "] y[" << p.y << "]}";
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

LonLat MapUnproject1(const LonLat& p0, const Vec& rp)
{
    double x_rad = rp.x / kEarthRadiusMean;
    double y_rad = rp.y / kEarthRadiusMean;
    double c = sqrt(x_rad * x_rad + y_rad * y_rad);
    double sin_c = sin(c);
    double cos_c = cos(c);
    double sin_lat0 = sin(p0.lat);
    double cos_lat0 = cos(p0.lat);

    double lat_rad;
    double lon_rad;
    if (fabs(c) > numeric_limits<double>::epsilon()) {
        lat_rad = asin(cos_c * sin_lat0 + (y_rad * sin_c * cos_lat0) / c);
        lon_rad = (p0.lon + atan2(x_rad * sin_c, c * cos_lat0 * cos_c - y_rad * sin_lat0 * sin_c));
    } else {
        lat_rad = p0.lat;
        lon_rad = p0.lon;
    }
    return { lon_rad, lat_rad };
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
    using Scalar = float;
    Scalar x_scale = cos(Scalar(p0.lat)) * Scalar(kEarthRadiusMean);
    Scalar y_scale = Scalar(kEarthRadiusMean);
    double x = Scalar(p1.lon - p0.lon) * x_scale;
    double y = Scalar(p1.lat - p0.lat) * y_scale;
    return { x, y };
}

LonLat MapUnproject4(const LonLat& p0, const Vec& rp)
{
    using Scalar = float;
    Scalar x_scale = cos(Scalar(p0.lat)) * Scalar(kEarthRadiusMean);
    Scalar y_scale = Scalar(kEarthRadiusMean);
    double lon = rp.x / x_scale + p0.lon;
    double lat = rp.y / y_scale + p0.lat;
    return { lon, lat };
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
        { 113.355309123, 23.152298456 },
        { 113.366920981229, 23.143200799852 },
        { 113.366920981223, 23.143200808156 }
    };
    for (size_t i = 0; i + 1 < pts.size(); ++i) {
        const auto& p1 = pts[i];
        const auto& p2 = pts[i + 1];
        auto p1_rad = ToRad(p1);
        auto p2_rad = ToRad(p2);
        auto pr1 = MapProject1(p1_rad, p2_rad);
        auto up1 = MapUnproject1(p1_rad, pr1);
        // auto pr2 = MapProject2(p1_rad, p2_rad);
        // auto pr3 = MapProject3(p1_rad, p2_rad);
        auto pr4 = MapProject4(p1_rad, p2_rad);
        auto up4 = MapUnproject4(p1_rad, pr4);
        auto& ref = pr1;
        cout << "p1: " << p1 << "\n"
             << "p2: " << p2 << "\n"
             << "pr1: " << pr1 << " " << Diff(ref, pr1) << "\n"
             //<< "pr2: " << pr2 << " " << Diff(ref, pr2) << "\n"
             //<< "pr3: " << pr3 << " " << Diff(ref, pr3) << "\n"
             << "pr4: " << pr4 << " " << Diff(ref, pr4) << "\n"
             << "up1:" << ToDeg(up1) << " " << up1 << "\n"
             << "up4:" << ToDeg(up4) << " " << up4 << "\n";
        cout << endl;
    }
}