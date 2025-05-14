template <typename T>
class SegmentTree {
public:
    enum class Type { MAX, MIN };

private:
    struct Node {
        T data, tag;
        bool hasTag;
        Node() : data(0), tag(0), hasTag(false) {}
    };

    size_t n;
    std::vector<Node> seg;
    Type        type;
    T           identity;
    std::function<T(T,T)> combine;

    T getVal(size_t id) const {
        return seg[id].hasTag ? seg[id].tag : seg[id].data;
    }

    void pull(size_t id) {
        seg[id].data = combine(getVal(id*2), getVal(id*2+1));
    }

    void push(size_t id) {
        if (!seg[id].hasTag) return;
        for (size_t c : {id*2, id*2+1}) {
            seg[c].tag    = seg[id].tag;
            seg[c].hasTag = true;
        }
        seg[id].data   = seg[id].tag;
        seg[id].hasTag = false;
    }

    T query(size_t ql, size_t qr, size_t l, size_t r, size_t id) {
        if (qr < l || ql > r) return identity;
        if (ql <= l && r <= qr) return getVal(id);
        push(id);
        size_t mid = (l + r) >> 1;
        return combine(
            query(ql, qr, l,     mid,   id*2),
            query(ql, qr, mid+1, r,     id*2+1)
        );
    }

    void update(T val, size_t ql, size_t qr, size_t l, size_t r, size_t id) {
        if (qr < l || ql > r) return;
        if (ql <= l && r <= qr) {
            seg[id].tag    = val;
            seg[id].hasTag = true;
            return;
        }
        push(id);
        size_t mid = (l + r) >> 1;
        update(val, ql, qr, l,     mid,   id*2);
        update(val, ql, qr, mid+1, r,     id*2+1);
        pull(id);
    }

public:
    SegmentTree() : n(0), type(Type::MAX), identity(0) {}

    void init(size_t nn, Type tp) {
        n = nn; 
        type = tp;
        if (type == Type::MAX) {
            identity = std::numeric_limits<T>::min();
            combine  = [](T a, T b){ return std::max(a,b); };
        } else {
            identity = std::numeric_limits<T>::max();
            combine  = [](T a, T b){ return std::min(a,b); };
        }
        seg.assign(n * 4, Node{});
    }

    T query(size_t ql, size_t qr) {
        return query(ql, qr, 0, n-1, 1);
    }
    void update(size_t ql, size_t qr, T val) {
        update(val, ql, qr, 0, n-1, 1);
    }
};
