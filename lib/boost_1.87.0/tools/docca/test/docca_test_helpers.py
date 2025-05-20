#
# Copyright (c) 2024 Dmitry Arkhipov (grisumbras@yandex.ru)
#
# Distributed under the Boost Software License, Version 1.0. (See accompanying
# file LICENSE_1_0.txt or copy at http://www.boost.org/LICENSE_1_0.txt)
#
# Official repository: https://github.com/boostorg/json
#

class MockXmlElem():
    def __init__(self, tag='compound'):
        self.tag = tag
        self.text = None
        self.tail = None

        self._attrs = dict()
        self._children = []

    def get(self, attr, default=None):
        return self._attrs.get(attr, default)

    def __len__(self):
        return len(self._children)

    def __getitem__(self, pos):
        return self._children[pos]

    def findall(self, tag):
        return [child for child in self._children if child.tag == tag]

    def find(self, tag):
        elems = self.findall(tag)
        if elems:
            return elems[0]

    def itertext(self):
        result = self.text or ''
        for child in self:
            result += child.itertext()
        result += self.tail or ''
        return result

    def append(self, child):
        return self._children.append(child)

    def clear(self):
        return self._children.clear()

    def update(self, **kw):
        return self._attrs.update(**kw)

def make_elem(contents):
    result = MockXmlElem(contents.get('tag'))

    prev = None
    for item in contents.get('items', []):
        if isinstance(item, str):
            assert not isinstance(prev, str)
            if prev is None:
                assert result.text is None
                result.text = item
            else:
                prev.tail = item
            prev = item
        else:
            child = make_elem(item)
            result.append(child)
            prev = child

    attrs = dict((
        (k,v) for (k,v) in contents.items()
        if k not in ('tag', 'items')
    ))
    result.update(**attrs)

    return result
