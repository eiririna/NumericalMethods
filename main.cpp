#include <bits/stdc++.h>
using namespace std;

int outputmatrix(vector<vector<double>> A) {
    int m = A.size(), n = A[0].size();
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++)
            cout << A[i][j] << " ";
        cout << endl;
    }
    cout << endl;
    return 0;
}

int outputvector(vector<double> a) {
    int n = a.size();
    for (int i = 0; i < n; i++)
        cout << a[i] << endl;
    cout << endl;
    return 0;
}

vector<vector<double>> matrixmult(vector<vector<double>> A, vector<vector<double>> B) {
    vector<vector<double>> AB;
    vector<double> ab;
    double c;
    for (int i = 0; i < A.size(); i++) {
        ab.clear();
        for (int j = 0; j < B[0].size(); j++) {
            c = 0;
            for (int k = 0; k < A[0].size(); k++)
                c += A[i][k] * B[k][j];
            ab.push_back(c);
        }
        AB.push_back(ab);
    }
    return AB;
}

vector<vector<double>> matrixtransp(vector<vector<double>> A) {
    vector<vector<double>> AT;
    vector<double> at;
    double h;
    for (int j = 0; j < A[0].size(); j++) {
        at.clear();
        for (int i = 0; i < A.size(); i++)
            at.push_back(A[i][j]);
        AT.push_back(at);
    }
    return AT;
}

vector<vector<double>> matrixsum(vector<vector<double>> A, vector<vector<double>> B) {
    for (int i = 0; i < A.size(); i++)
        for (int j = 0; j < A[0].size(); j++)
            A[i][j] += B[i][j];
    return A;
}

vector<vector<double>> matrixmultconst(vector<vector<double>> A, double c) {
    for (int i = 0; i < A.size(); i++)
        for (int j = 0; j < A[0].size(); j++)
            A[i][j] *= c;
    return A;
}

double vectornorma(vector<vector<double>> A) {
    double s = 0;
    for (int i = 0; i < A.size(); i++)
        s += pow(A[i][0], 2);
    return sqrt(s);
}

vector<double> gaussmethod(vector<vector<double>> Ab) {
    vector<vector<double>> A, L, U, P;
    vector<double> a, l, b, x, y, Pb;
    int i, j, n = Ab.size();
    double h;

    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++)
            a.push_back(Ab[i][j]);
        A.push_back(a);
        a.clear();
    }

    for (i = 0; i < n; i++)
        b.push_back(Ab[i][n]);

    a.clear();
    for (i = 0; i < n; i++) {
        for (j = 0; j < n; j++)
            if (i == j)
                a.push_back(1);
            else
                a.push_back(0);
        P.push_back(a);
        a.clear();
    }

    U = A;
    a.clear();
    l.clear();
    for (int k = 0; k < n - 1; k++) {
        i = 0;
        if (U[k][k] == 0) {
            i++;
            while (U[k + i][k] == 0)
                i++;
            for (j = 0; j < n; j++) {
                a.push_back(U[k][j]);
                U[k][j] = U[k + i][j];
                U[k + i][j] = a[j];
            }
            a.clear();
            for (j = k; j < n; j++) {
                a.push_back(P[k][j]);
                P[k][j] = P[k + i][j];
                P[k + i][j] = a[j - k];
            }
            a.clear();
        }
        for (i = 0; i < k; i++)
            l.push_back(0);
        l.push_back(1);
        for (i = k + 1; i < n; i++) {
            h = U[i][k] / U[k][k];
            for (j = k; j < n; j++)
                U[i][j] -= h * U[k][j];
            l.push_back(h);
        }
        L.push_back(l);
        l.clear();
    }
    for (i = 0; i < n - 1; i++)
        l.push_back(0);
    l.push_back(1);
    L.push_back(l);

    for (i = 0; i < n - 1; i++)
        for (j = i + 1; j < n; j++) {
            h = L[i][j];
            L[i][j] = L[j][i];
            L[j][i] = h;
        }

    for (i = 0; i < b.size(); i++) {
        Pb.push_back(0);
        for (j = 0; j < b.size(); j++)
            Pb[i] += P[i][j] * b[j];
    }

    for (i = 0; i < n; i++) {
        y.push_back(Pb[i]);
        for (j = 0; j < i; j++)
            y[i] -= L[i][j] * y[j];
    }

    for (i = n - 1; i > -1; i--) {
        x.push_back(y[i]);
        for (j = i + 1; j < n; j++)
            x[n - 1 - i] -= U[i][j] * x[n - 1 - j];
        x[n - 1 - i] /= U[i][i];
    }
    a.clear();
    for (i = 0; i < x.size(); i++)
        a.push_back(x[x.size() - 1 - i]);
    x = a;

    return x;
}

vector<vector<double>> zeidelmethod(vector<vector<double>> Ab) {
    vector<vector<double>> A, b, C, d, x;
    vector<double> a;
    int n = Ab.size(), i;
    double a0, e = 1e-6;

    for (int i = 0; i < n; i++) {
        a.clear();
        for (int j = 0; j < n; j++)
            a.push_back(Ab[i][j]);
        A.push_back(a);
    }
    for (int i = 0; i < n; i++)
        b.push_back({ Ab[i][n] });

    b = matrixmult(matrixtransp(A), b);
    A = matrixmult(matrixtransp(A), A);

    for (int k = 0; k < n - 1; k++) {
        i = 0;
        if (A[k][k] == 0) {
            i++;
            a.clear();
            while (A[k + i][k] == 0)
                i++;
            for (int j = 0; j < n; j++) {
                a.push_back(A[k][j]);
                A[k][j] = A[k + i][j];
                A[k + i][j] = a[j];
            }
            a0 = b[k][0];
            b[k][0] = b[k + i][0];
            b[k + i][0] = a0;
        }
    }

    for (int i = 0; i < n; i++) {
        a.clear();
        for (int j = 0; j < n; j++)
            if (i == j)
                a.push_back(0);
            else
                a.push_back(-A[i][j] / A[i][i]);
        C.push_back(a);
        d.push_back({ b[i][0] / A[i][i] });
    }

    x = d;
    while (vectornorma(matrixsum(matrixmult(A, x), matrixmultconst(b, -1))) > e)
        for (int i = 0; i < n; i++) {
            x[i][0] = d[i][0];
            for (int j = 0; j < n; j++)
                x[i][0] += C[i][j] * x[j][0];
        }

    return x;
}

int main() {
    vector<vector<double>> Ab;
    vector<double> a;
    string s, s0;
    int pos;

    do {
        getline(cin, s);
        s.append(" ");
        a.clear();
        do {
            pos = s.find(" ");
            s0.append(s, 0, pos);
            s.erase(0, pos + 1);

            const char* c = s0.c_str();
            a.push_back(strtod(c, NULL));
            s0.clear();
        } while (s != "");
        Ab.push_back(a);
    } while (Ab.size() < a.size() - 1);
    cout << endl;

    outputvector(gaussmethod(Ab));
    outputmatrix(zeidelmethod(Ab));

    return 0;
}