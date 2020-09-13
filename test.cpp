#include <iostream>
#include <vector>
#include <list>
#include <algorithm>
#include <queue>
#include <unordered_map>

using namespace std;

long merge(vector<int> &arr, int si, int mid, int ei)
{
    vector<int> res;
    int i = si, j = mid + 1;
    int n = mid + 1, m = ei + 1;
    int count = 0;

    while (i < n && j < m)
    {
        if (arr[i] < arr[j])
            res.push_back(arr[i++]);
        else
        {
            res.push_back(arr[j++]);
            count += n - i;
        }
    }

    while (i < n)
        res.push_back(arr[i++]);
    while (j < m)
        res.push_back(arr[j++]);

    for (int k = si; k < ei + 1; k++)
        arr[k] = res[k - si];

    return count;
}
int mergeSort(vector<int> &arr, int si, int ei)
{
    if (si == ei)
        return 0;

    long mid = (si + ei) / 2;
    int lcount = mergeSort(arr, si, mid);
    int rcount = mergeSort(arr, mid + 1, ei);

    int myCount = merge(arr, si, mid, ei);

    return lcount + rcount + myCount;
}
long countInversions(vector<int> &arr)
{
    return mergeSort(arr, 0, arr.size() - 1);
}

void solve()
{
    int t;
    cin >> t;
    while (t-- > 0)
    {
        int n;
        cin >> n;
        vector<int> arr(n);
        for (int i = 0; i < n; i++)
            cin >> arr[i];
        cout << countInversions(arr) << endl;
    }
}
int main()
{
    // g++ test.cpp -o out && out < input.txt > output.txt
    solve();
    return 0;
}