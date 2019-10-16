## Adding A New Doxygen Page

1. Have a tag for your .md's first title header such as `{\#my_md}`
2. Take your .md and place it in an appropriate folder (or make one)
3. In `pretty_stuff/layout.xml`, in the `navindex`:
4. Add a new tab within an existing tab
    - For instance, to add a new user guide, under the tab with title "Contributor Guide"
      add a new tab such as
```
<tab type="usergroup" visible="yes" title="My Cool New Guide" url="@ref my_md" intro=""/>
```
