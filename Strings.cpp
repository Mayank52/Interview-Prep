#include <iostream>
#include <vector>
#include <list>
#include <unordered_map>
#include <map>
#include <algorithm>

using namespace std;

// Group words with same set of characters
vector<vector<string>> groupWords(vector<string> &words)
{
    map<vector<bool>, vector<string>> mp; //char freq: words

    //make a map of (char present:words)
    for (int i = 0; i < words.size(); i++)
    {
        vector<bool> freq(26, 0);
        for (char c : words[i])
            freq[c - 'a'] = true;

        mp[freq].push_back(words[i]);
    }

    vector<vector<string>> res;

    //push all words with same same character set into the result
    for (auto ele : mp)
        res.push_back(ele.second);

    for (vector<string> &word : res)
    {
        for (string str : word)
            cout << str << " ";
        cout << endl;
    }

    return res;
}

// 1347. Minimum Number of Steps to Make Two Strings Anagram
int minSteps(string s, string t)
{
    vector<int> freq(26);
    //count the  freq of characters in str1
    for (int i = 0; i < s.size(); i++)
        freq[s[i] - 'a']++;

    int count = 0;
    //check for the characters in str2, if they have same freq
    for (int i = 0; i < t.size(); i++)
    {
        if (freq[t[i] - 'a'] == 0)
            count++;
        else
            freq[t[i] - 'a']--;
    }

    return count;
}

// 890. Find and Replace Pattern
bool match(string &word, string &pattern)
{
    unordered_map<char, char> mp;

    //1.check that all characters of word are mapped to unique values
    //map all values in word to pattern
    for (int i = 0; i < word.size(); i++)
    {
        //if it present in map but mapped to some other value then return false
        if (mp.find(word[i]) != mp.end())
        {
            if (mp[word[i]] != pattern[i])
                return false;
        }
        //else add it to the map
        mp[word[i]] = pattern[i];
    }

    //2. check if all values are mapped to unique keys
    //check if more than one value has been mapped to same value i.e. more than one key has same value
    vector<bool> seen(26, false);
    for (auto ele : mp)
    {
        //if the value has already been mapped, then return false
        if (seen[ele.second - 'a'])
            return false;
        //seeing it for first time
        seen[ele.second - 'a'] = true;
    }

    return true;
}
vector<string> findAndReplacePattern(vector<string> &words, string pattern)
{
    vector<string> res;
    for (string &word : words)
    {
        if (match(word, pattern))
            res.push_back(word);
    }

    return res;
}

// Generate all binary strings from given pattern
void binaryStrings(string &str, int idx, string res)
{
    if (idx == str.size())
    {
        cout << res << endl;
        return;
    }

    if (str[idx] == '?')
    {
        binaryStrings(str, idx + 1, res + '0');
        binaryStrings(str, idx + 1, res + '1');
    }
    else
    {
        binaryStrings(str, idx + 1, res + str[idx]);
    }
}

// 76. Minimum Window Substring
string minWindow(string s, string t)
{
    int n = s.size();
    int m = t.size();

    if (m > n)
        return "";

    unordered_map<char, int> mp1;
    unordered_map<char, int> mp2;

    //make freq map of T string
    for (char c : t)
        mp1[c]++;

    int i = -1, j = -1, currCount = 0, reqCount = t.size();
    int si = 0, ei = 0, minLen = s.size();

    while (true)
    {
        bool f1 = false;
        bool f2 = false;

        //acquire
        while (i < n - 1 && currCount < reqCount)
        {
            f1 = true;

            i++;
            mp2[s[i]]++;

            if (mp2[s[i]] <= mp1[s[i]])
                currCount++;
        }

        //remove
        while (j < i && currCount == reqCount)
        {
            f2 = true;

            j++;
            mp2[s[j]]--;

            if (mp2[s[j]] < mp1[s[j]])
                currCount--;
        }

        //update minLen
        if (i - j < minLen)
        {
            si = j;
            ei = i;
            minLen = i - j;
        }

        if (!f1 && !f2)
            break;
    }

    if (j == -1)
        return "";

    return s.substr(si, ei - si + 1);
}

// 49. Group Anagrams
//Approach 1-Group by sorted strings
//Approach 2- Make a Hashmap of {freqMap: strings}
vector<vector<string>> groupAnagrams(vector<string> &strs)
{
    vector<vector<string>> res;
    map<vector<int>, vector<string>> mp; //freq map: strings

    //group all strings with same freq map in the map
    for (string &str : strs)
    {
        //make freq map
        vector<int> freq(26, 0);
        for (char ch : str)
        {
            freq[ch - 'a']++;
        }
        //add the string to the map
        mp[freq].push_back(str);
    }

    for (auto ele : mp)
        res.push_back(ele.second);

    return res;
}

// 438. Find All Anagrams in a String
//Approach - Compare both freq maps, for each window, Time: O(n*26)= O(n)
vector<int> findAnagrams(string s, string p)
{
    int n = s.size();
    int m = p.size();

    if (n < m)
        return {};

    vector<int> res;
    vector<int> smp(26), pmp(26);

    //make freq for p string
    for (char ch : p)
        pmp[ch - 'a']++;

    int i;

    //freq map of first window of s
    for (i = 0; i < m; i++)
        smp[s[i] - 'a']++;

    if (smp == pmp)
        res.push_back(0);

    //other windows
    for (; i < n; i++)
    {
        smp[s[i] - 'a']++;
        smp[s[i - m] - 'a']--;
        //compare the current window with p map
        if (pmp == smp)
            res.push_back(i - m + 1);
    }

    return res;
}

// Check if two strings are k-anagrams or not
bool areKAnagrams(string str1, string str2, int k)
{
    if (str1.size() != str2.size())
        return false;

    int freq[26] = {0};
    int count = 0;

    //make freq map of str1
    for (int i = 0; i < str1.size(); i++)
        freq[str1[i] - 'a']++;

    for (int i = 0; i < str2.size(); i++)
    {
        //if it is in freq map, then decrease it
        if (freq[str2[i] - 'a'] > 0)
            freq[str2[i] - 'a']--;
        //else increase the count
        else
            count++;
    }

    return count > k ? false : true;
}

// Check if binary representations of two numbers are anagram
/*
Just count the number of set bits in both numbers
if they are equal then they are anagrams
*/
int setBitCount(int n)
{
    //Using Kernighan’s Algorithm
    int count = 0;
    while (n)
    {
        n &= (n - 1);
        count++;
    }
    return count;
}
bool checkBinaryAnagram(int a, int b)
{
    if (setBitCount(a) == setBitCount(b))
        return true;
    return false;
}

// Make largest palindrome by changing at most K-digits
string maximumPalinUsingKChanges(string str, int k)
{
    string palin = str;

    // Iinitialize l and r by leftmost and
    // rightmost ends
    int l = 0;
    int r = str.length() - 1;

    //  first try to make string palindrome
    while (l < r)
    {
        // Replace left and right character by
        // maximum of both
        if (str[l] != str[r])
        {
            palin[l] = palin[r] = max(str[l], str[r]);
            k--;
        }
        l++;
        r--;
    }

    // If k is negative then we can't make
    // string palindrome
    if (k < 0)
        return "Not possible";

    l = 0;
    r = str.length() - 1;

    while (l <= r)
    {
        // At mid character, if K>0 then change
        // it to 9
        if (l == r)
        {
            if (k > 0)
                palin[l] = '9';
        }

        // If character at lth (same as rth) is
        // less than 9
        if (palin[l] < '9')
        {
            /* If none of them is changed in the 
               previous loop then subtract 2 from K 
               and convert both to 9 */
            if (k >= 2 && palin[l] == str[l] &&
                palin[r] == str[r])
            {
                k -= 2;
                palin[l] = palin[r] = '9';
            }

            /*  If one of them is changed in the previous 
                loop then subtract 1 from K (1 more is 
                subtracted already) and make them 9  */
            else if (k >= 1 && (palin[l] != str[l] ||
                                palin[r] != str[r]))
            {
                k--;
                palin[l] = palin[r] = '9';
            }
        }
        l++;
        r--;
    }

    return palin;
}

// Lexicographically first palindromic string
/*
Make Freq map of the string
for palindrome the freq of at most 1 char can be odd, if not then return not possible
iterate through the freq map,and put each char at first and last pos, until its count becomes 0
if it is odd then save it for the middle
*/
string firstPalindromicString(string &str)
{
    string res = str;
    int freq[26] = {0};
    for (int i = 0; i < str.size(); i++)
        freq[str[i] - 'a']++;

    int midEle = -1, l = 0, r = str.size() - 1;
    for (int i = 0; i < str.size(); i++)
    {
        //if current freq is odd, and we also have a prev off freq, then not possible
        if (freq[i] % 2 != 0 && midEle != -1)
            return "Not Possible";

        else
        {
            //if even freq, add it to current left and right index, and update the freq
            if (freq[i] % 2 == 0)
            {
                while (freq[i] > 0)
                {
                    char ch = i + 'a';
                    res[l] = res[r] = ch;
                    l++;
                    r--;
                    freq[i] -= 2;
                }
            }
            //if it is odd
            else
            {
                //set current ele as the mid ele
                midEle = i;
                while (freq[i] > 1)
                {
                    char ch = i + 'a';
                    res[l] = res[r] = ch;
                    l++;
                    r--;
                    freq[i] -= 2;
                }
            }
        }
    }

    //set the middle element if length was odd
    if (midEle != -1)
        res[str.size() / 2] = (char)(midEle + 'a');

    return res;
}

// Longest Non-palindromic substring
/*
Check if the complete string is palindrome or not
If it is, then remove a  character from the start or end of string gives the longest non palindrome
If it is not a palindrome, then the whole string is the ans

Only exception is if all characters of string are same, then there is no non palindromic substring
*/
string nonPalindromeSubstring(string &str)
{
    int n = str.size();
    int ch = str[0], f1 = 1, f2 = 1; //f1: to check if all chars in string are same, f2 to check if palindrome
    for (int i = 0; i < n / 2; i++)
    {
        if (str[i] != ch)
            f1 = 0;
        if (str[i] != str[n - i - 1])
            f2 = 0;
    }

    //all same
    if (f1 == 1)
        return "Not Possible";
    //complete palindrome
    if (f2 == 1)
        return str.substr(1, str.size());
    //not palindrome
    else
        return str;
}



void solve()
{
}
int main()
{
    solve();
    return 0;
}