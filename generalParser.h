//
//  Created by Federico Pagnozzi on 28/11/14.
//  Copyright (c) 2014 Federico Pagnozzi. All rights reserved.
//  This file is distributed under the BSD 2-Clause License. See LICENSE.TXT
//  for details.
#ifndef GENERALPARSER_H
#define GENERALPARSER_H
#include "emilibase.h"

#include <exception>
#include <string>
#include <vector>
#include <map>
#include <set>
#include <sstream>
#include <functional>
#include <iomanip>

namespace prs
{
void emili_header();
void info();
void check(const char *t, const char* message);
void printTab(const char* string);
void printTab(const std::string& string);
void incrementTabLevel();
void decrementTabLevel();

struct TabLevel {
    TabLevel() { incrementTabLevel(); }
    ~TabLevel() { decrementTabLevel(); }
};

struct ParsingError : std::exception {
    std::string msg;
    ParsingError() {}
    ParsingError(std::string msg_) : msg(msg_) {}

    const char* what() const throw() override {
        return msg.c_str();
    }
};

class TokenManager
{
protected:
    const char** tokens;
    int numberOfTokens;
    int currentToken;
public:
    TokenManager(const char** tokens, int numberOfTokens):tokens(tokens),numberOfTokens(numberOfTokens),currentToken(0) { }

    bool hasNext() const;
    void next();
    void operator ++(int); // same as next
    const char* nextToken(); // skip one token and return it, or null if no token left

    const char* peek() const; // return current token or " "
    const char* operator *() const; // same as peek

    /**
     * @example 50 -> 50
     * @example bla -> Error
     * @example -> Error
     */
    int getInteger() throw(ParsingError);

    /**
     * @example 50 -> 50.0
     * @example 5.2 -> 5.2
     * @example 2e3 -> 2000.0
     * @example bla -> Error
     * @example -> Error
     */
    float getDecimal() throw(ParsingError);

    /**
     * @example (N=20) percent 50 -> 10
     * @example (N=20) % 50 -> 10
     * @example (N=20) %50 -> 10
     * @example (N=20) 50% -> 10
     * @example (N=20) 50 -> 50
     * @example (N=20) -> Error
     * @example (N=20) bla -> Error
     * @example (N=20) % -> Error
     * @example (N=20) % bla -> Error
     */
    float getPercent(float N) throw(ParsingError);

    bool checkToken(const std::string &token);
    bool checkToken(const char*);

    const char* tokenAt(int i) const; // get i'th token, or null
    bool peekIs(const char*) const;
    bool peekIs(const std::string &) const;

    /**
     * @brief if(the next token is an integer) { res = the integer; return true } else { return false }
     * @param res the location to write the integer
     * @return true if an integer was set to res
     */
    bool checkInteger(int &res);

    /**
     * @brief if(the next token is a decimal) { res = the decimal; return true } else { return false }
     * @param res the location to write the decimal
     * @return true if an decimal was set to res
     */
    bool checkDecimal(float &res);

    /**
     * @brief if(the next token is a decimal) { res = the decimal; return true } else { return false }
     * @param res the location to write the decimal
     * @return true if an decimal was set to res
     */
    bool checkDecimal(double &res);
};

struct ErrorExpected : ParsingError {
    ErrorExpected(prs::TokenManager& tm, std::string name, std::vector<std::string> const& tokensRaw) {
        auto startsWith = [](std::string s, std::string g) {
            return s.substr(0, g.size()) == g;
        };

        auto endsWith = [](std::string s, std::string g) {
            return s.substr(s.size() - g.size(), g.size()) == g;
        };

        auto getSlice = [](std::string s, int i, int j) -> std::string {
            if(i < 0) i = s.size() + i;
            if(j < 0) j = s.size() + j;
            return s.substr(i, j-i);
        };

        bool everyLine = endsWith(tm.peek(), "?");
        std::string prefix = getSlice(tm.peek(), 0, -1);

        std::ostringstream oss;
        oss << "'" << tm.peek() << "' -> ERROR a " << name << " is expected : ";

        std::vector<std::string> tokensFiltered;
        std::vector<std::string> const * tokens = &tokensRaw;

        if(everyLine) {
            for(auto x : tokensRaw)
                if(startsWith(x, prefix))
                    tokensFiltered.push_back(x);
            tokens = &tokensFiltered;
        }

        if(tokens->size() == 0) {
            if(everyLine)
                oss << "<" << '\n' << ">";
            else
                oss << "<>";
        } else {
            if(! everyLine) {
                auto it = tokens->begin();
                oss << "<" << *it++;
                while(it != tokens->end())
                    oss << " | " << *it++;
                oss << ">";
            } else {
                for(auto x : *tokens)
                    oss << "\n" << (prefix.empty() ? "" : "... ") << x.substr(prefix.size());
            }
        }
        msg = oss.str();
    }
};

/**
 * @brief ArgParser parse command line args into a dict(str: (int|float|bool))
 *
 * examples:
 *
 *  // numbers with defaults
 *  addInt("a", 10);
 *  addInt("b", 20);
 *  addFloat("x", 1.5);
 *  (1) a 10
 *  (2) { a 10 }
 *  (3) a 10 c ...
 *  -> {a:10, b:20, x:1.5}
 *  (4) { a 10 c 20 }
 *  -> ERROR "c" is not a valid param
 *
 *  in (3), parsing stops at "c" because it is not a valid param
 *  in (4), error because "c" is not a valid param
 *
 *  // bools using {param} or no-{param}
 *  addBool("premium", false);
 *  addBool("valid", true);
 *  (1)
 *  -> {premium: false, valid: true}
 *  (2) valid premium
 *  -> {premium: true, valid: true}
 *  (3) no-valid premium
 *  -> {premium: true, valid: false}
 *
 *  // required
 *  addIntRequired("a", 5);
 *  (1)
 *  (2) { }
 *  -> ERROR "a" required
 */
struct ArgParser {
    std::map<std::string, int> dataInt;
    std::map<std::string, bool> dataBool;
    std::map<std::string, float> dataFloat;
    std::map<std::string, std::string> dataString;
    typedef std::function<void(TokenManager&)> Func;
    std::map<std::string, Func> dataFunc;
    std::vector<std::string> order;
    std::set<std::string> given;
    std::vector<std::string> required;
    std::string name;

    enum Tag {
        REQUIRED
    };

    ArgParser() {
        name = "param";
    }

    ArgParser(std::string parent) {
        name = parent + " param";
    }

    // add

    void addIntRequired(std::string s) {
        addInt(s, 0);
        required.push_back(s);
    }

    void addStringRequired(std::string s) {
        addString(s, "");
        required.push_back(s);
    }

    void addBoolRequired(std::string s) {
        addBool(s, false);
        required.push_back(s);
    }

    void addFloatRequired(std::string s) {
        addFloat(s, 0.0f);
        required.push_back(s);
    }

    void addFuncRequired(std::string s, Func d) {
        addFunc(s, d);
        required.push_back(s);
    }

    bool exists(std::string s) {
        return dataInt.count(s) || dataBool.count(s) || dataFloat.count(s) || dataString.count(s) || dataFunc.count(s);
    }

    void addInt(std::string s, int d) {
        if(exists(s))
            throw std::invalid_argument(" argument " + s + " already exist !");
        dataInt[s] = d;
        order.push_back(s);
    }

    void addString(std::string s, std::string d) {
        if(exists(s))
            throw std::invalid_argument(" argument " + s + " already exist !");
        dataString[s] = d;
        order.push_back(s);
    }

    void addFloat(std::string s, float d) {
        if(exists(s))
            throw std::invalid_argument(" argument " + s + " already exist !");
        dataFloat[s] = d;
        order.push_back(s);
    }

    void addBool(std::string s, bool d) {
        if(exists(s))
            throw std::invalid_argument(" argument " + s + " already exist !");
        dataBool[s] = d;
        order.push_back(s);
    }

    void addFunc(std::string s, Func d) {
        if(exists(s))
            throw std::invalid_argument(" argument " + s + " already exist !");
        dataFunc[s] = d;
        order.push_back(s);
    }

    // synonyms for add

    void addInt(std::string s, Tag t) {
        if(t == REQUIRED)
            addIntRequired(s);
    }

    void addBool(std::string s, Tag t) {
        if(t == REQUIRED)
            addBoolRequired(s);
    }

    void addFloat(std::string s, Tag t) {
        if(t == REQUIRED)
            addFloatRequired(s);
    }

    void addString(std::string s, Tag t) {
        if(t == REQUIRED)
            addStringRequired(s);
    }

    void addFunc(std::string s, Func d, Tag t) {
        if(t == REQUIRED)
            addFuncRequired(s, d);
    }

    // get

    int Int(std::string s) {
        if(! dataInt.count(s))
            throw std::invalid_argument("Missing int " + s);
        return dataInt[s];
    }

    float Float(std::string s) {
        if(! dataFloat.count(s))
            throw std::invalid_argument("Missing int " + s);
        return dataFloat[s];
    }

    bool Bool(std::string s) {
        if(! dataBool.count(s))
            throw std::invalid_argument("Missing bool " + s);
        return dataBool[s];
    }

    std::string const& String(std::string s) {
        if(! dataString.count(s))
            throw std::invalid_argument("Missing bool " + s);
        return dataString[s];
    }

    bool has(std::string s) {
        if(! exists(s))
            throw std::invalid_argument("Argument " + s + " does not exist");
        return given.count(s);
    }

    // methods

    void operator()(prs::TokenManager& tm) {
        parse(tm);
    }

    /**
     * @brief parseSequence so that multiple { } are like a same { }
     * @example { a 5 } { b 6 } will be the same as { a 5 b 6 }
     * Caution: must begin with "{" to be parsed
     * @example a 5 b 6 will not be parsed
     */
    void parseSequence(prs::TokenManager& tm){
        while (tm.peekIs("{"))
            parse(tm, false);
        testRequired();
    }

    void parse(prs::TokenManager& tm, bool callTestRequired = true) {
        /*
         * Simple Algo :
         * for(;;)
         *     if(tm.checkToken("A"))
         *         A = tm.getInteger()
         *     else if()...
         *     else if(stop) break
         *     else if(!peekIs("}"));
         * More over, may start/end with {}
         */
        // Imagine data* attributes have "a" and "b" but no "x"
        // a 5 b 2 => OK
        // a 5 x 1 b 2 => OK but b is not read...
        // { a 5 b 2 } => OK
        // { a 5 b 2 x 1 } => NOT OK ("x" found expected "}")
        // { a 5 x 1 b 2 } => NOT OK ("x" found expected "}")

        auto endsWith = [](std::string s, std::string g) {
            return s.substr(s.size() - g.size(), g.size()) == g;
        };

        bool hasBrace = tm.checkToken("{");

        for(;;) {
            bool found = false;

            if(!found)
            for(auto& p : dataInt) {
                if(tm.checkToken(p.first) || tm.checkToken(p.first + ":")) {
                    tm.checkToken(":");
                    p.second = tm.getInteger();
                    found = true;
                    given.insert(p.first);
                    break;
                }
            }

            if(! found)
            for(auto& p : dataString) {
                if(tm.checkToken(p.first) || tm.checkToken(p.first + ":")) {
                    tm.checkToken(":");
                    auto n = tm.nextToken();
                    if(!n)
                        throw std::invalid_argument("expected string token");
                    p.second = n;
                    found = true;
                    given.insert(p.first);
                    break;
                }
            }

            if(!found)
            for(auto& p : dataFloat) {
                if(tm.checkToken(p.first) || tm.checkToken(p.first + ":")) {
                    tm.checkToken(":");
                    p.second = tm.getDecimal();
                    found = true;
                    given.insert(p.first);
                    break;
                }
            }

            if(!found)
            for(auto& p : dataBool) {
                if(tm.checkToken(p.first)) {
                    p.second = true;
                    found = true;
                    given.insert(p.first);
                    break;
                }
                else if(tm.checkToken("no-" + p.first)) {
                    p.second = false;
                    found = true;
                    given.insert(p.first);
                    break;
                }
            }

            if(!found)
            for(auto& p : dataFunc) {
                if(tm.checkToken(p.first) || tm.checkToken(p.first + ":")) {
                    tm.checkToken(":");
                    p.second(tm);
                    found = true;
                    given.insert(p.first);
                    break;
                }
            }

            if(! found) {
                if(endsWith(tm.peek(), "?"))
                    throw ErrorExpected(tm, name, order);
                else if(!hasBrace)
                    break;
                else if(tm.peekIs("}"))
                    break;
                else
                    throw ErrorExpected(tm, name, order);
            }
        }

        // assert(given <= required);
        if(callTestRequired)
            testRequired();

        if(hasBrace && ! tm.checkToken("}"))
            throw ErrorExpected(tm, "}", {"}"});
    }

    /**
     * @brief throw an error if any required parameter is not given
     */
    void testRequired() {
        std::vector<std::string> r;
        for(auto s : required)
            if(! given.count(s))
                r.push_back(s);

        if(! r.empty()) {
            if(r.size() == 1)
                throw ParsingError("Argument '" + r.front() + "' is required");
            else {
                auto it = r.begin();
                std::ostringstream oss;
                oss << *it++;
                while(it != r.end())
                    oss << ", " << *it++;

                throw ParsingError("Arguments '" + oss.str() + "' are required");
            }
        }
    }

    /**
     * @brief Sort by type, then name
     */
    void printTypeAlpha() {
        TabLevel l;
        int m = 0;
        for(auto x : order)
            m = std::max(m, (int)x.size());

        for(auto& p : dataInt) {
            std::ostringstream oss;
            oss << "|" << std::setw(m) << p.first << ": " << p.second << std::endl;
            printTab(oss.str());
        }

        for(auto& p : dataString) {
            std::ostringstream oss;
            oss << "|" << std::setw(m) << p.first << ": " << p.second << std::endl;
            printTab(oss.str());
        }

        for(auto& p : dataFloat) {
            std::ostringstream oss;
            oss << "|" << std::setw(m) << p.first << ": " << p.second << std::endl;
            printTab(oss.str());
        }

        for(auto& p : dataBool) {
            std::ostringstream oss;
            oss << "|" << std::setw(m) << p.first << ": " << std::boolalpha << p.second << std::endl;
            printTab(oss.str());
        }

        // no func print (already done when called)
    }

    void printInline(std::string msg) {
        std::ostringstream oss;
        oss << msg << "(";
        const char* f = "";
        for(auto x : order) {
            oss << f;
            f = ", ";

            if(dataInt.count(x))
                oss << x << " = " << dataInt[x];
            else if(dataFloat.count(x))
                oss << x << " = " << dataFloat[x];
            else if(dataBool.count(x))
                oss << x << " = " << std::boolalpha << dataBool[x];
            else if(dataString.count(x))
                oss << x << " = " << dataString[x];

            // no func print (already done when called)
        }

        oss << ")";
        printTab(oss.str());
    }

    void print(std::string msg) {
        printTab(msg);
        print();
    }

    /**
     * @brief Insertion order
     */
    void print() {
        TabLevel l;
        int m = 0;
        for(auto x : order)
            m = std::max(m, (int)x.size());
        for(auto x : order) {
            std::ostringstream oss;

            if(dataInt.count(x))
                oss << std::setw(m) << x << ": " << dataInt[x];
            else if(dataFloat.count(x))
                oss << std::setw(m) << x << ": " << dataFloat[x];
            else if(dataBool.count(x))
                oss << std::setw(m) << x << ": " << std::boolalpha << dataBool[x];
            else if(dataString.count(x))
                oss << std::setw(m) << x << ": " << dataString[x];
            else
                continue; // no func print (already done when called)

            printTab(oss.str());
        }
    }

};

struct CheckBrac {
    bool has;
    TokenManager& tm;

    CheckBrac(TokenManager& tm_) : tm(tm_) {
        has = tm.checkToken("[");
    }

    ~CheckBrac() {
        if(has && ! tm.checkToken("]"))
            throw ErrorExpected(tm, "]", {"]"});
    }
};

class AlgoBuilder
{
public:
    virtual std::string availableProblems() const{ return std::string("Iamabstract!");}
    virtual std::vector<std::string> availableProblemsList() const { return {}; }
public:
    virtual bool isParsable(std::string& problem)=0 ; // should be const std::string &
    virtual emili::LocalSearch* buildAlgo(prs::TokenManager& tm) {return nullptr;}
    virtual std::string info() {return std::string("Iamabstract!");}
    virtual bool operator ==(const AlgoBuilder& b);

    /**
     * may return null
     */
    virtual emili::Problem* getInstance() {
        return nullptr;
    }
};

class NoSearch : public std::exception {
    const char* what() const throw() {
        return "NoSearch";
    }
};

class GeneralParser
{
protected:
    std::vector< AlgoBuilder* > builders;
    TokenManager tm;
public:        
    GeneralParser(const char** tokens, int numberOfTokens):tm(tokens,numberOfTokens) { }
    GeneralParser(TokenManager tokenmanager):tm(tokenmanager) { }
    virtual emili::LocalSearch* parseParams();
    virtual void registerBuilder(AlgoBuilder* builder);
    virtual void removeBuilder(AlgoBuilder* builder);
    bool noSearch = false;
};
}

#endif
