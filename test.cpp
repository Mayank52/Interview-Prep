#include <iostream>
#include <vector>
#include <list>
#include <algorithm>
#include <queue>
#include <unordered_map>

using namespace std;

int maxProfit(int shares, int profit, int K, int idx, vector<int> &prices, vector<vector<vector<int>>> &dp)
{

    if (K == 0 || idx == prices.size())
    {
        return dp[idx][K][shares] = profit;
    }

    if (dp[idx][K][shares] != -1)
        return dp[idx][K][shares];

    int maxp = 0;
    //buy
    if (shares == 0)
    {
        maxp = max(maxp, maxProfit(1, profit - prices[idx], K, idx + 1, prices, dp));
    }

    //sell
    else
    {
        maxp = max(maxp, maxProfit(0, profit + prices[idx], K - 1, idx + 1, prices, dp));
    }

    //do nothing
    return dp[idx][K][shares] = max(maxp, maxProfit(shares, profit, K, idx + 1, prices, dp));

}
int maxProfit_01(int k, vector<int> &prices)
{
    int n = prices.size();
    vector<vector<vector<int>>> dp(n + 1, vector<vector<int>>(k + 1, vector<int>(2, -1)));

    int ans = maxProfit(0, 0, k, 0, prices, dp);

    for (int i = 0; i < n + 1; i++)
    {
        for (int j = 0; j < k + 1; j++)
        {
            cout << "(" << dp[i][j][0] << "," << dp[i][j][1] << ")"
                 << " ";
        }
        cout << endl;
    }

    return ans;
}

void solve()
{
    vector<int> prices{7, 2, 3, 4, 5};
    maxProfit_01(2, prices);
}
int main()
{
    // g++ Arrays.cpp -o out && ./out < input.txt > output.txt
    solve();
    return 0;
}