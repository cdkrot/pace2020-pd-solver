// Dmitry _kun_ Sayutin (2019)

#ifdef DEBUG_MODE
#define _GLIBCXX_DEBUG
#define LOCAL
#endif

#include <bits/stdc++.h>
#include <execinfo.h>
#include <unistd.h>

using std::cin;
using std::cout;
using std::cerr;

using std::vector;
using std::map;
using std::unordered_map;
using std::array;
using std::set;
using std::string;
using std::queue;

//#define map unordered_map

using std::pair;
using std::make_pair;

using std::tuple;
using std::make_tuple;
using std::get;

using std::min;
using std::abs;
using std::max;
using std::swap;

using std::shared_ptr;
using std::unique_ptr;

using std::unique;
using std::sort;
using std::generate;
using std::reverse;
using std::min_element;
using std::max_element;

#define SZ(vec)         int((vec).size())
#define ALL(data)       data.begin(),data.end()
#define RALL(data)      data.rbegin(),data.rend()
#define TYPEMAX(type)   std::numeric_limits<type>::max()
#define TYPEMIN(type)   std::numeric_limits<type>::min()

#define ensure(cond) if (not (cond)) {fprintf(stderr, "Failed: %s", #cond); exit(1);}

#undef assert

void print_trace() { // adapted from gnu docs
    void* array[10];
    size_t size;
    char** strings;
    size_t i;

    size = backtrace(array, 10);
    strings = backtrace_symbols(array, size);

    for (i = 0; i < size; i++)
        cerr << strings[i] << "\n";

    free(strings);
}

void attach_debugger() {
#ifdef DEBUG_MODE
    volatile bool is_debugged = false;
    while (not is_debugged) {
        usleep(1);
    }
#endif
}

#ifdef DEBUG_MODE
#define assert(X) {if (not (X)) {cerr << "assertion failed\n"; print_trace(); while(1) {}; exit(1);}}
#else
#define assert(X) {}
#endif

extern "C" void abort_handler(int signal_number) {
#ifdef DEBUG_MODE
    assert(false);
#else
#endif
}


struct instance_t {
    int n;
    int m;
    vector<vector<int>> graph;

    int large() const {
        return 5 * n;
    }
};

struct result_t {
    int depth;
    int root_id;
    vector<int> tree_parents;
};

instance_t parse(std::istream& in) {
    string s;
    int n = -1, m = -1;
    
    vector<vector<int>> adj;
    
    while (std::getline(in, s)) {
        if (s.empty() or s[0] == 'c')
            continue;

        if (s[0] == 'p') {
            ensure(s.substr(0, 6) == "p tdp ");
            s = s.substr(6, 100000);

            const char* ptr = s.data();
            n = strtoll(ptr, const_cast<char**>(&ptr), 10);
            m = strtoll(ptr, const_cast<char**>(&ptr), 10);
            
            adj.resize(n);
        } else {
            const char* ptr = s.data();
            
            int a = strtoll(ptr, const_cast<char**>(&ptr), 10) - 1;
            int b = strtoll(ptr, const_cast<char**>(&ptr), 10) - 1;
            adj[a].push_back(b);
            adj[b].push_back(a);
        }
    }

    return instance_t {n, m, adj};
}

void print_result(std::ostream& os, instance_t& instance, result_t& result) {
    os << result.depth << "\n";
    for (int i = 0; i < instance.n; ++i)
        os << (i == result.root_id ? 0 : result.tree_parents[i] + 1) << "\n";
}

template <typename T>
using static_set = vector<T>;

template <typename T, typename Container>
static_set<T> make_set(const Container& cont) {
    vector<T> res;
    for (const T& elem: cont)
        res.push_back(elem);

    sort(ALL(res));
    assert(unique(ALL(res)) == res.end());
    return res;
}

template <typename T>
static_set<T> remove(static_set<T> st, T& x) {
    auto it = std::lower_bound(ALL(st), x);
    if (it != st.end())
        st.erase(it);

    return st;
}


template <typename T>
static_set<T> add(static_set<T> st, T& x) {
    auto it = std::lower_bound(ALL(st), x);
    assert(it == st.end() or *it != x);
    st.insert(it, x);

    return st;
}

template <typename T>
static_set<T> set_diff(const static_set<T>& a, const static_set<T>& b) {
    static_set<T> res;
    std::set_difference(ALL(a), ALL(b), std::back_inserter(res));
    return res;
}

template <typename T>
int id_in(const static_set<T>& st, T object) {
    auto rs = std::lower_bound(ALL(st), object);

    if (rs == st.end() or *rs != object)
        return -1;
    return rs - st.begin();
}

enum color_t {
    BLACK = 1,
    WHITE = 0,
    GROUNDSET_COLOR = 2,
};

struct marked_forest_t {
    static_set<int> vertices;

    map<int, char> labels;
    map<int, vector<int>> children;
    map<int, int> parent; // set to -1 for roots


    int find_empty_id(int minimum = 0) const {
        int ans = minimum;
        while (id_in(vertices, ans) != -1)
            ++ans;

        return ans;
    }

    void replace_vertex(int x, int y) {
        if (x == y)
            return;
        
        vertices = add(remove(vertices, x), y);

        labels[y] = get_label(x);
        labels.erase(labels.find(x));

        children[y] = std::move(children[x]);
        children.erase(children.find(x));

        parent[y] = std::move(parent[x]);
        parent.erase(parent.find(x));

        for (auto& elem: parent)
            if (elem.second == x)
                elem.second = y;

        for (auto& elem: children)
            for (int& ch: elem.second)
                if (ch == x)
                    ch = y;
    }

    void add_vertex(int v, int par, color_t color) {
        vertices = add(vertices, v);
        
        labels[v] = color;

        children[v] = {};
        parent[v] = par;

        if (par != -1)
            children[par].push_back(v);
    }
    
    template <typename Func>
    void for_each_root(Func func) const {
        for (int v: vertices)
            if (get_parent(v) == -1)
                func(v);
    }

    int get_label(int v) const {
        auto it = labels.find(v);
        assert(it != labels.end());
        return it->second;
    }

    int get_parent(int v) const {
        auto it = parent.find(v);
        assert(it != parent.end());
        return it->second;
    }

    const vector<int>& get_children(int v) const {
        auto it = children.find(v);
        assert(it != children.end());
        return it->second;
    }
    
    int get_depth(int v) const {
        int ans = 0;
        while (get_parent(v) != -1)
            v = get_parent(v), ++ans;
        return ans;
    }

    int level_ancestor(int v, int k) const {
        for (int i = 0; i < k; ++i)
            v = get_parent(v);
        return v;
    }
    
    bool is_in_closure(int v, int u) const {
        int dv = get_depth(v);
        int du = get_depth(u);

        if (dv > du)
            v = level_ancestor(v, dv - du);
        if (du > dv)
            u = level_ancestor(u, du - dv);

        return (v == u);
    }
    

};

template <typename Key, typename Val>
struct trie_map {
    trie_map() {
    }

    template <typename Cont>
    std::optional<Val>& get_if_exists(const Cont& cont) {
        node_t* cur = root.get();

        for (Key key: cont) {
            if (not cur->go[key].get())
                return nullval;

            cur = cur->go[key].get();
        }

        return cur->val;
    }
    
    template <typename Cont>
    std::optional<Val>& operator[](const Cont& cont) {
        node_t* cur = root.get();

        for (Key key: cont) {
            if (not cur->go[key].get())
                cur->go[key] = std::move(std::make_unique<node_t>());

            cur = cur->go[key].get();
        }

        return cur->val;
    }

    template <typename Cont>
    Val& get_or_create(const Cont& cont) {
        auto& ref = (*this)[cont];

        if (ref == std::nullopt)
            ref = Val();

        return *ref;
    }
    
private:
    struct node_t {
        std::optional<Val> val;
        map<Key, unique_ptr<node_t>> go;
    };

    std::optional<Val> nullval = std::nullopt;
    unique_ptr<node_t> root = std::make_unique<node_t>();
};

struct partial_decomp_t;
std::optional<std::map<int, int>> establish_mapping(const partial_decomp_t& a, const partial_decomp_t& b);

enum class parenting_type {
    RESTRICT,
    JOIN,
    EXTEND,
    NONE
};

struct decomp_parental_info {
    shared_ptr<const partial_decomp_t> parent;
    shared_ptr<const partial_decomp_t> parent2;
    parenting_type type = parenting_type::NONE;
};

struct partial_decomp_t: public std::enable_shared_from_this<partial_decomp_t> {
    partial_decomp_t() = default;
    partial_decomp_t(const marked_forest_t& forest, const static_set<int>& x_set,
                     const decomp_parental_info& info): forest(forest), x_set(x_set), info(info) {}
    
    marked_forest_t forest;
    static_set<int> x_set;

    decomp_parental_info info;
    
    bool operator==(const partial_decomp_t& other) const {
        auto res = establish_mapping(*this, other);
        if (res == std::nullopt)
            return false;

        for (int v: forest.vertices)
            if (forest.get_label(v) != other.forest.get_label((*res)[v]))
                return false;

        return true;
    }

    shared_ptr<const partial_decomp_t> clone() const {
        partial_decomp_t res = *this;
        return std::make_shared<const partial_decomp_t>(std::move(res));
    }

    partial_decomp_t clone_to_value() const {
        partial_decomp_t res = *this;
        return res;
    }

    partial_decomp_t& operator=(const partial_decomp_t& other) = default;

    partial_decomp_t(const partial_decomp_t& other) = default;
};

struct isomorphism_checker_t {
    vector<int> add_forest(const marked_forest_t& forest, bool account_marks) {
        std::function<int(int)> dfs = [&](int v) {
            vector<int> ch_ids;
            for (int ch: forest.children.find(v)->second)
                ch_ids.push_back(dfs(ch));

            if (account_marks)
                ch_ids.push_back(-100 + forest.get_label(v));

            return id_of(ch_ids);
        };

        vector<int> lst;
        for (int v: forest.vertices)
            if (forest.get_parent(v) == -1)
                lst.push_back(dfs(v));

        sort(ALL(lst));
        return lst;
    }

    vector<int> add_pd_noting_groundset(const partial_decomp_t& decomp, bool add_new=true) {
        vector<int> res(SZ(decomp.x_set), -1);
        
        std::function<int(int)> dfs = [&](int v) {
            vector<int> ch_ids;
            for (int ch: decomp.forest.children.find(v)->second)
                ch_ids.push_back(dfs(ch));

            if (decomp.forest.get_label(v) == GROUNDSET_COLOR)
                ch_ids.push_back(-100 - v);

            int new_id = (add_new ? id_of(ch_ids) : id_of_ifexists(ch_ids));
            if (new_id == -1)
                throw -1;

            if (decomp.forest.get_label(v) == GROUNDSET_COLOR)
                res[id_in(decomp.x_set, v)] = new_id;

            return new_id;
        };

        try {
            decomp.forest.for_each_root(dfs);
        } catch (int ex) {
            return vector<int> {-322}; // nonexistent ID
        };

        return res;
    }
    
private:
    int id_of(vector<int>& ids) {
        sort(ALL(ids));
        auto& res = id[ids];
        if (not res.has_value())
            res = next_id++;

        return res.value();
    }

    int id_of_ifexists(vector<int>& ids) {
        sort(ALL(ids));

        auto val = id.get_if_exists(ids);
        if (val == std::nullopt)
            return -1;
        return val.value();
    }

    
    trie_map<int, int> id;
    int next_id = 0;
};

std::optional<std::map<int, int>> establish_mapping(const partial_decomp_t& a, const partial_decomp_t& b) {
    if (a.x_set != b.x_set or SZ(a.forest.vertices) != SZ(b.forest.vertices))
        return std::nullopt;
    
    std::map<int, int> mapping;
    set<int> used;
    for (int v: a.x_set) {
        int v1 = v;
        int v2 = v;

        while (v1 != -1 and v2 != -1) {
            if (mapping.count(v1) and mapping[v1] != v2)
                return std::nullopt;

            if (mapping.count(v1)) {
                v1 = v2 = -1;
                break;
            }

            if (used.count(v2))
                return std::nullopt;
            
            mapping[v1] = v2;
            used.insert(v2);
            
            v1 = a.forest.get_parent(v1);
            v2 = b.forest.get_parent(v2);
        }

        if (v1 != -1 or v2 != -1)
            return std::nullopt; // even depth is different
    }

    return mapping;
}

struct decomp_set {
    vector<shared_ptr<const partial_decomp_t>> lst;
    isomorphism_checker_t checker;
    set<vector<int>> graphset;
    
    auto begin() {
        return lst.begin();
    }
    
    auto end() {
        return lst.end();
    }

    auto begin() const {
        return lst.cbegin();
    }
    
    auto end() const {
        return lst.cend();
    }

    bool empty() const {
        return lst.empty();
    }
    
    void add(const partial_decomp_t& a) {
        for (int v: a.forest.vertices)
            assert(a.forest.get_parent(v) != v);
        
        // for (const auto& x: lst)
        //     if (*x == a)
        //         return;

        auto desc = checker.add_forest(a.forest, true);
        if (graphset.count(desc))
            return;
        
        lst.push_back(a.shared_from_this());
    }

    void add(const decomp_set& other) {
        for (auto elem: other)
            add(*elem);
    }
    
    template <typename Func>
    void sequence_to(Func func, decomp_set& res) {
        for (auto elem: *this)
            func(*elem, res);
    }
    
    template <typename Func>
    decomp_set sequence(Func func) {
        decomp_set res;
        
        sequence_to(func, res);

        return res;
    }
};

shared_ptr<partial_decomp_t> restrict(const partial_decomp_t& decomp,
                                                 static_set<int> xprime,
                                                 static_set<int> y,
                                                 static_set<int> z) {

    partial_decomp_t result;
    result.x_set = xprime;

    for (int v: xprime)
        assert(id_in(decomp.x_set, v) != -1);
    
    std::function<bool(int)> dfs = [&](int v) {
        vector<int> alive_ch;

        for (int go: decomp.forest.get_children(v))
            if (dfs(go))
                alive_ch.push_back(go);

        if (alive_ch.empty() and id_in(xprime, v) == -1)
            return false;

        color_t mycolor;
        if (id_in(xprime, v) != -1)
            mycolor = GROUNDSET_COLOR;
        else if (id_in(y, v) != -1)
            mycolor = BLACK;
        else if (id_in(z, v) != -1)
            mycolor = WHITE;
        else
            mycolor = (color_t)decomp.forest.get_label(v);
        
        result.forest.add_vertex(v, -1, mycolor);

        for (int u: alive_ch) {
            result.forest.parent[u] = v;
            result.forest.children[v].push_back(u);
        }
        
        return true;
    };

    decomp.forest.for_each_root(dfs);
    
    result.info.parent = decomp.shared_from_this();
    result.info.type = parenting_type::RESTRICT;
    return std::make_shared<partial_decomp_t>(std::move(result));
}

void extend(const partial_decomp_t& a_, int vert, decomp_set& out,
            int height_limit, const instance_t& inst) {
    
    partial_decomp_t a = a_.clone_to_value();
    
    if (id_in(a.forest.vertices, vert) != -1) {
        assert(false);
        
        a.forest.replace_vertex(vert, a.forest.find_empty_id(inst.large()));
        assert(id_in(a.x_set, vert) == -1);
    }

    static_set<int> new_x_set = add(a.x_set, vert);

    auto add_if_success = [&](marked_forest_t& changed_forest) {
        for (int u: inst.graph[vert])
            if (id_in(new_x_set, u) != -1 and not changed_forest.is_in_closure(vert, u))
                return;
        
        decomp_parental_info info;
        info.parent = a_.shared_from_this();
        info.type = parenting_type::EXTEND;
        
        shared_ptr<partial_decomp_t> res = std::make_shared<partial_decomp_t>
            (changed_forest, new_x_set, info);
        out.add(*res.get());
    };
    
    for (int v: a.forest.vertices) {
        if (a.forest.get_label(v) == WHITE) {
            auto changed_forest = a.forest;
            changed_forest.replace_vertex(v, vert);

            changed_forest.labels[vert] = GROUNDSET_COLOR;

            add_if_success(changed_forest);
        }

        if (true or a.forest.get_label(v) == GROUNDSET_COLOR) {
            int depth_cur = a.forest.get_depth(v);

            for (int i = 1; i + depth_cur < height_limit; ++i) {
                auto changed_forest = a.forest;

                int last = v;
                for (int j = 1; j <= i; ++j) {
                    int newv = changed_forest.find_empty_id(inst.large());
                    changed_forest.add_vertex(newv, last, WHITE);

                    last = newv;
                }

                changed_forest.replace_vertex(last, vert);
                changed_forest.labels[vert] = GROUNDSET_COLOR;

                add_if_success(changed_forest);
            }
        }
    }

    // with "-1" parent
    if (1) {
        for (int i = 0; i < height_limit; ++i) {
            auto changed_forest = a.forest;
            
            int last = -1;
            for (int j = 0; j <= i; ++j) {
                int newv = changed_forest.find_empty_id(inst.large());
                changed_forest.add_vertex(newv, last, WHITE);
                
                last = newv;
            }
            
            changed_forest.replace_vertex(last, vert);
            changed_forest.labels[vert] = GROUNDSET_COLOR;
            
            add_if_success(changed_forest);
        }   
    }
}

shared_ptr<const partial_decomp_t> join(const partial_decomp_t& a, const partial_decomp_t& b) {
    auto mapping = establish_mapping(a, b);

    if (mapping == std::nullopt)
        return nullptr;

    partial_decomp_t res = a.clone_to_value();
    for (int v: a.forest.vertices) {
        if (a.forest.get_label(v) == GROUNDSET_COLOR)
            continue;

        int cnt = int(a.forest.get_label(v) == BLACK) + int(b.forest.get_label((*mapping)[v]) == BLACK);
        if (cnt == 2)
            return nullptr;

        res.forest.labels[v] = (cnt == 1 ? BLACK : WHITE);
    }

    res.info.parent = a.shared_from_this();
    res.info.parent2 = b.shared_from_this();
    res.info.type = parenting_type::JOIN;
    
    return std::make_shared<partial_decomp_t>(std::move(res));
}


// assume that tree is rooted at v=0, and all edges go from smaller numbered end to larger numbered end.
std::optional<result_t> solve(vector<static_set<int>> bags, vector<vector<int>> tw_go, instance_t& instance, int height_limit) {
    vector<decomp_set> results(SZ(bags));

    shared_ptr<partial_decomp_t> trivial_ptr = std::make_shared<partial_decomp_t>();
    const partial_decomp_t& trivial = *trivial_ptr;
    
    for (int i = SZ(bags) - 1; i >= 0; --i) {
        auto adapt = [&](const partial_decomp_t& elem_, decomp_set& res) {
            auto elem = elem_.clone();
            
            static_set<int> to_remove, to_add;
            
            for (int v: elem->x_set)
                if (id_in(bags[i], v) == -1)
                    to_remove.push_back(v);

            for (int v: bags[i])
                if (id_in(elem->x_set, v) == -1)
                    to_add.push_back(v);

            
            if (not to_remove.empty())
                elem = restrict(*elem, set_diff(elem->x_set, to_remove), to_remove, static_set<int> {});

            assert(SZ(elem->forest.vertices) == SZ(elem->forest.parent));
            
            decomp_set cur;
            cur.add(*elem);
            
            for (int a: to_add) {
                cur = cur.sequence([&](const partial_decomp_t& b, decomp_set& rs) {
                    extend(b, a, rs, height_limit, instance);
                });
            }

            res.add(cur);
        };
        
        decomp_set& cur = results[i];

        if (tw_go[i].empty()) {
            cur.add(trivial);

            cur = cur.sequence(adapt);
        } else {
            cur = std::move(results[tw_go[i][0]]).sequence(adapt);

            for (int j = 1; j < SZ(tw_go[i]); ++j) {
                auto tmp = results[tw_go[i][j]].sequence(adapt);

                // consider only pairs with bijection
                
                isomorphism_checker_t iso;

                if (not (SZ(cur.lst) < SZ(tmp.lst)))
                    swap(cur, tmp);

                trie_map<int, vector<int>> mp;
                decomp_set res;
                for (int i = 0; i < SZ(cur.lst); ++i) {
                    mp.get_or_create(iso.add_pd_noting_groundset(*cur.lst[i])).push_back(i);
                }
                
                for (auto& a: tmp) {
                    auto& sublst = mp.get_if_exists(iso.add_pd_noting_groundset(*a, false));
                    if (sublst == std::nullopt)
                        continue;
                    
                    for (int b: *sublst) {
                        auto joinres = join(*a, *cur.lst[b]);

                        if (joinres)
                            res.add(*joinres);
                        
                        if (i == 0 and j + 1 == SZ(tw_go[i]) and not res.empty())
                            break;
                    }

                    if (i == 0 and j + 1 == SZ(tw_go[i]) and not res.empty())
                        break;
                }

                cur = std::move(res);
            }
        }
        
        cerr << height_limit << "," << i << ": " << results[i].lst.size() << "\n";

        if (results[i].empty())
            return std::nullopt;
    }
    
    if (results[0].empty())
        return std::nullopt;

    int n = instance.n;
    vector<int> parents(n, -1);

    const partial_decomp_t& decomp = *results[0].lst[0];

    marked_forest_t forest;

    std::map<int, int> true_location_of;
    std::map<int, int> by_location_val;
    
    std::function<void(const partial_decomp_t&, std::map<int, int>&)> recover = [&](const partial_decomp_t& decomp, std::map<int, int>& mapping) {
        std::function<void(int)> work = [&](int v) {
            if (not mapping.count(v)) {
                int parent = decomp.forest.get_parent(v);
                int mapped_parent = (parent == -1 ? -1 : mapping[parent]);

                int mapped_id;
                if (id_in(forest.vertices, v) == -1) {
                    mapped_id = v;
                } else {
                    mapped_id = forest.find_empty_id(instance.large());
                }
                
                forest.add_vertex(mapped_id, mapped_parent, BLACK);
                mapping[v] = mapped_id;
            }

            for (int u: decomp.forest.get_children(v))
                work(u);
        };

        decomp.forest.for_each_root(work);

        for (int v: decomp.forest.vertices)
            if (v < instance.large()) {
                true_location_of[v] = mapping[v];
                by_location_val[mapping[v]] = v;
            }
        
        if (decomp.info.type == parenting_type::NONE)
            return;

        if (decomp.info.type == parenting_type::EXTEND) {
            const partial_decomp_t& parent_decomp = *decomp.info.parent;
            
            for (auto it = mapping.begin(); it != mapping.end();)
                if (id_in(parent_decomp.forest.vertices, it->first) == -1) {
                    mapping.erase(it++);
                } else {
                    ++it;
                }
            
            recover(parent_decomp, mapping);
            return;
        }

        if (decomp.info.type == parenting_type::RESTRICT) {
            recover(*decomp.info.parent, mapping);
            return;
        }

        if (decomp.info.type == parenting_type::JOIN) {
            auto submap = *establish_mapping(*decomp.info.parent2, *decomp.info.parent);

            for (auto& elem: submap)
                elem.second = mapping[elem.second];

            recover(*decomp.info.parent, mapping);            
            recover(*decomp.info.parent2, submap);
            
            return;
        }
    };

    std::map<int, int> tmpmap;
    recover(decomp, tmpmap);
    
    for (int v = 0; v < n; ++v) {
        int u = forest.get_parent(true_location_of[v]);
        while (u != -1 and not by_location_val.count(u))
            u = forest.get_parent(u);

        parents[v] = (u == -1 ? -1 : by_location_val[u]);
    }
    
    assert(std::count(ALL(parents), -1) == 1);
    int root_id = std::find(ALL(parents), -1) - parents.begin();
    return result_t {height_limit, root_id, parents};
}

vector<vector<int>> parse_list_of_lists(int* data) {
    int len = data[0];
    assert(data[1] == -2);

    data += 2;
    
    vector<vector<int>> res;

    for (int i = 0; i < len; ++i) {
        res.push_back(vector<int> {});
        
        while (*data != -1)
            res.back().push_back(*data++);
        data++;
    }

    return res;
}

extern "C" {
    void python_enter_point_ng(int n, int m, int* edges, int* buckets, int* tree) {
        signal(SIGABRT, abort_handler);
        
        instance_t instance;
        instance.n = n;
        instance.m = m;

        instance.graph.resize(n);
    
        for (int i = 0; i < m; ++i) {
            int v = edges[i] / n;
            int u = edges[i] % n;

            instance.graph[v].push_back(u);
            instance.graph[u].push_back(v);
        }

        auto bags = parse_list_of_lists(buckets);
        auto treedecomp = parse_list_of_lists(tree);

        for (auto& bag: bags)
            sort(ALL(bag));
        
        cerr << "BAGS\n";
        for (int i = 0; i < SZ(bags); ++i) {
            cerr << i << ": ";
            for (auto ch: bags[i])
                cerr << ch << ",";
            cerr << "\n";
        }

        cerr << "Decomp\n";
        for (int i = 0; i < SZ(treedecomp); ++i) {
            cerr << i << ": ";
            for (auto ch: treedecomp[i])
                cerr << ch << ",";
            cerr << "\n";
        }

        cerr << "Solving:\n";

        int H;
        std::optional<result_t> res;
        for (H = 1; H <= 20; ++H) {
            if (res = solve(bags, treedecomp, instance, H))
                break;
        }
        if (not res)
            assert(false);
        
        cerr << "H was " << H << "\n";
        print_result(cout, instance, *res);
        exit(0);
    }
}

// #ifdef AS_MAIN
// int main() {
//     std::iostream::sync_with_stdio(false);
//     cin.tie(nullptr);
//     cout.tie(nullptr);

//     // code here
//     instance_t instance = parse(cin);
        
//     auto result = solve(instance);
//     print_result(cout, instance, result);
    
//     return 0;
// }
// #endif
