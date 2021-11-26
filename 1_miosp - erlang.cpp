#include <iostream>
#include <fstream>
#include <string>
#include <vector>


#include "pbPlots.hpp"
#include "supportLib.hpp"


struct FourVal
{
    double w, x, y, z;
    

    FourVal(double w, double x,
            double y = 0., double z = 0.) {
        this->w = w;
        this->x = x;
        this->y = y;
        this->z = z;
    }      
};

struct FourValText {
    std::string wText, xText, yText, zText;

    FourValText(std::string wText, std::string xText,
        std::string yText = "", std::string zText = "") {
        this->wText = wText;
        this->xText = xText;
        this->yText = yText;
        this->zText = zText;
    }
};

//z1


void save2csv(std::string fileName, std::vector<FourVal> values, FourValText col) {
    std::ofstream dataOut;

    dataOut.open(fileName + ".csv", std::ios_base::app);
    if (dataOut.is_open()) {
        dataOut << col.wText << ";" << col.xText << col.yText << col.zText << std::endl;
        for (std::vector<FourVal>::iterator it = values.begin(); it != values.end(); ++it) {
            //std::cout << it->first << std::setw(10) << it->second << std::endl;
            dataOut << it->w << ";" << it->x << ";" << it->y << ";" << it->z << std::endl;
        }
        dataOut.close();
        std::cout << "save to csv done!\n";
    }
    else std::cout << "save error!\n";
    return;
}

long long int factorial(int n)
{
    if (n == 0) return 1;
    else return n * factorial(n - 1);
}

std::vector<std::pair<double, double>> probability(double a_min, double a_max, double a_step, double V) {
    std::vector<std::pair<double, double>> values;
    std::vector<FourVal> a;
    for (double i = a_min; i <= a_max; i += a_step) {
        double sum = 0;
        double A = i * V;
        double up = pow(A, V) / factorial(V);


        for (int j = 0; j <= V; j++) {
            double down = (pow(A, j) / factorial(j));
            sum += down;

        }
        values.push_back(std::make_pair(A, up / sum));
        a.push_back(FourVal(A, up / sum));

    }
    save2csv("erlangFile", a, FourValText("A", "val"));
    return values;
}

long double elRec(double A, double V) {
    if (V == 0) return 1;
    else {
        double con = A * elRec(A, V - 1);
        return con / (V + con);
    }
}


std::vector<std::pair<double, double>> probability2(double a_min, double a_max, double a_step, double V) {
    std::vector<std::pair<double, double>> values; 
    std::vector<FourVal> a;
    for (double i = a_min; i <= a_max; i += a_step) {
        double sum = 0;
        double A = i * V;

        values.push_back(std::make_pair(A, elRec(A, V)));
        a.push_back(FourVal(A, elRec(A, V)));
    }
    save2csv("erlangFile2", a,FourValText("A", "val"));
    return values;
}



void drawPNG(std::vector<std::pair<double, double>> val, std::string name) {
    std::vector<double> xs; //= { -2, -1, 0, 1, 2 };
    std::vector<double> ys;// = { 2, -1, -2, -1, 2 };

    for (std::vector<std::pair<double, double>>::iterator it = val.begin(); it != val.end(); ++it) {
        xs.push_back(it->first);
        ys.push_back(it->second);
    }

    RGBABitmapImageReference* imageReference = CreateRGBABitmapImageReference();
    DrawScatterPlot(imageReference, 600, 400, &xs, &ys);

    std::vector<double>* pngData = ConvertToPNG(imageReference->image);
    WriteToFile(pngData, "plot-" + name + ".png");
    DeleteImage(imageReference->image);
}

std::vector<std::pair<double, double>> predictA(double blocking, double V, double max_error = 0.0005, double astep = 0.1) {  // dla max_error = 0.005 znajduje takie same wartoœci jak na erlang.com
    std::vector<std::pair<double, double>> values;
    std::vector<FourVal> a;
    double A = 0, prob_Val = 0, abs_error = 0;
    
    do {
        A += astep;
        prob_Val = elRec(A, V);
        abs_error = abs(prob_Val - blocking);
        values.push_back(std::make_pair(A, abs_error));
        a.push_back(FourVal(A, abs_error));
        std::cout << A << " " << abs_error << std::endl;
    } while (abs_error >= max_error);
    std::cout << "A wynosi: " << A << "   " << std::endl;
    save2csv("predictA", a,FourValText("A", "abs_error"));
    return values;
}


std::vector<std::pair<double, double>> predictV(double blocking, double A, double max_error = 0.01, double v_step = 1) {  // dla max_error = 0.01 znajduje takie same wartoœci jak na erlang.com
    std::vector<std::pair<double, double>> values;
    std::vector<FourVal> a;
    double V = 0, prob_Val = 0, abs_error = 0;
    do {
        V += v_step;
        prob_Val = elRec(A, V);
        abs_error = abs(prob_Val - blocking);
        values.push_back(std::make_pair(V, abs_error));
        a.push_back(FourVal(V, abs_error));
        std::cout << V << " " << abs_error << std::endl;
    } while (abs_error >= max_error);
    std::cout << "V wynosi: " << V << "   " << std::endl;
    save2csv("predictV", a, FourValText("V", "abs_error"));
    return values;
}

std::vector<std::pair<double, double>> palm_jacobaeus(double V, double x = 1, double a_min = 1, double a_max = 10, double a_step = 0.5) {
    std::vector<std::pair<double, double>> values;
    std::vector<FourVal> a;
    double A, erl, H_x;
    while (a_min <= a_max) {
        A = a_min * V;
        erl = elRec(A, V);
        H_x = erl / elRec(A, V - x);

        std::cout << a_min << " : " << erl << " : " << H_x << std::endl;
        a.push_back(FourVal(a_min, erl, H_x));
        a_min += a_step;
        
    }
    save2csv("palm_jacobaeus", a, FourValText("V", "abs_error", "H_x"));
    return values;
}

int main()
{
    //drawPNG(probability(0.2, 1.3, 0.1, 1), "pic2");
    //drawPNG(probability2(0.2, 1.3, 0.1, 1), "pic3");

    predictA(0.5, 30);
    predictV(0.5, 30); 
    
    palm_jacobaeus(15);
    return 0;
}




