"""
Minimal Git Interface

This module implements a minimal git interface. We only use it to find
out the SHA1 of tree objects of the current branch.

Limitations:

* Readonly!

* No support for packed anything. You can be sure that ``'HEAD'``
  works, but other branches might not.
"""


import os
import re
import zlib

from .logger import logger
from .git_ignore import is_ignored


BLOB_RE = re.compile(r'(?P<type>[a-z]*) (?P<size>[0-9]*)\0(?P<content>.*)', flags=re.DOTALL)

BLOB_COMMIT_TREE_RE = re.compile(r'tree (?P<sha1>[a-f0-9]{40,40})')


def Blob(filename):
    """
    Load the git object with given sha1 hash
    """
    with open(filename, 'rb') as f:
        blob = zlib.decompress(f.read())
    match = BLOB_RE.match(blob)
    if not match:
        raise ValueError('file is not a git object')
    if blob.startswith('commit'):
        return BlobCommit(match)
    if blob.startswith('tree'):
        return BlobTree(match)
    if blob.startswith('blob'):
        return BlobFile(match)
    raise ValueError('unsupported blob: ' + repr(blob))


class BlobABC(object):
    
    def __init__(self, match):
        self._type = match.group('type')
        self._size = int(match.group('size'))
        self._content = match.group('content')

    def __repr__(self):
        s = [self.__class__.__name__ + ':']
        for line in self._content.splitlines():
            if len(line.strip()) == 0:
                continue
            s.append('    ' + line.rstrip())
        return '\n'.join(s)


class BlobCommit(BlobABC):
    
    @property
    def tree(self):
        """
        Return the tree sha1 that is being committed
        
        EXAMPLES::

            >>> commit = git.get_symbolic_ref('HEAD')
            >>> type(commit)
            <class 'sage_pkg.tiny_git.BlobCommit'>
            >>> commit.tree    # doctest: +SKIP
            'f1795efaddd86273a715baf26a80c908a132a035'
            >>> len(commit.tree) == 40 and int(commit.tree, 16) >= 0
            True
        """
        match = BLOB_COMMIT_TREE_RE.match(self._content)
        return match.group('sha1')

    
class BlobTree(BlobABC):

    MODE_NORMAL = '100644'
    MODE_EXEC = '100755'
    MODE_DIR = '40000'

    def ls(self):
        """
        Iterate over the tree content
        
        EXAMPLES::

            >>> commit = git.get_symbolic_ref('HEAD')
            >>> tree = git.get(commit.tree)
            >>> type(tree)
            <class 'sage_pkg.tiny_git.BlobTree'>
            >>> for item in tree.ls():   # doctest: +ELLIPSIS
            ...     print(item)
            ('100644', '.gitignore', ...
        """
        pos = 0
        while pos < len(self._content):
            pos_space = self._content.find(' ', pos + 1)
            mode = self._content[pos:pos_space]
            pos_zero = self._content.find('\0', pos_space + 1)
            name = self._content[pos_space + 1:pos_zero]
            sha1_binary = self._content[pos_zero + 1:pos_zero + 21]
            sha1 = ''.join('{0:02x}'.format(ord(c)) for c in sha1_binary)
            yield (mode, name, sha1)
            pos = pos_zero + 21

    def __repr__(self):
        """
        Iterate over the tree content
        
        Override the generic repr since the tree object contains binary sha1's.

        EXAMPLES::

            >>> commit = git.get_symbolic_ref('HEAD')
            >>> git.get(commit.tree)   # doctest: +ELLIPSIS
            BlobTree:
                100644 ...
        """
        s = [self.__class__.__name__ + ':']
        for mode, name, sha1 in self.ls():
            s.append('    {0:>6} {1} {2}'.format(mode, sha1, name))
        return '\n'.join(s)

    def ls_dirs(self):
        """
        Iterater over the sub-trees
        
        EXAMPLES::

            >>> commit = git.get_symbolic_ref('HEAD')
            >>> tree = git.get(commit.tree)
            >>> type(tree)
            <class 'sage_pkg.tiny_git.BlobTree'>
            >>> for item in tree.ls_dirs():   # doctest: +ELLIPSIS
            ...     print(item)
            ('40000', 'build', ...
        """
        for dirent in self.ls():
            if dirent[0] == self.MODE_DIR:
                yield dirent

    def get(self, filename):
        """
        Return the sha1 of a file in the tree

        EXAMPLES::

            >>> commit = git.get_symbolic_ref('HEAD')
            >>> tree = git.get(commit.tree)
            >>> sha1 = tree.get('.gitignore')
            >>> sha1  # doctests: +SKIP
            'a9b6be08742e31b728f04f6f89c4b93f28ac4b92'
            >>> len(sha1) == 40 and int(sha1, 16) >= 0
            True
        """
        for mode, name, sha1 in self.ls():
            if name == filename:
                return sha1
        raise ValueError('file is not in the tree')

    def _compare_files(self, dirname, verbose=False):
        """
        Helper to just compare the files
        """


class BlobFile(BlobABC):
    pass
    



class GitRepository(object):

    def __init__(self, dot_git):
        if not os.path.isdir(dot_git):
            raise ValueError('directory does not exist')
        self.dot_git = dot_git

    def __repr__(self):
        """
        Return a string representation

        EXAMPLES::

            >>> git    # doctest: +ELLIPSIS
            Git repository at .../.git
        """
        return 'Git repository at {0}'.format(self.dot_git)

    def get(self, sha1):
        """
        Load the blob with given sha1

        INPUT:

        - ``sha1`` -- 40-digit hex number as string. The sha1 of the
          object to load. Must be an unpacked blob.
        """
        logger.debug('Loading blob with sha1 = %s', sha1)
        filename = os.path.join(self.dot_git, 'objects', sha1[0:2], sha1[2:])
        return Blob(filename)

    def get_symbolic_ref(self, symbolic_ref):
        """
        Load the blob with given symbolic ref

        INPUT:

        - ``symbolic_ref`` -- string. A symbolic ref like ``'HEAD'``,
          ``'FETCH_HEAD'``, ....

        EXAMPLES::

            >>> git.get_symbolic_ref('HEAD')     # doctest: +ELLIPSIS
            BlobCommit:
                tree ...
        """
        logger.debug('Looking up symbolic ref = %s', symbolic_ref)
        filename = os.path.join(self.dot_git, symbolic_ref)
        with open(filename, 'r') as f:
            ref = f.read()
        before, after = ref.split('ref: ')
        if len(before) > 0:
            raise ValueError('not a symbolic ref')
        return self.get_branch(after.rstrip())

    def get_branch(self, branch):
        """
        Load the blob that the branch head points to

        INPUT:

        - ``branch`` -- a full branch name, starting with ``'refs/'``.
        """
        logger.debug('Looking up branch = %s', branch)
        if not branch.startswith('refs/'):
            raise ValueError('branch name must be full')
        filename = os.path.join(*([self.dot_git] + branch.split(os.path.sep)))
        if os.path.exists(filename):
            # Ref is unpacked
            with open(filename, 'r') as f:
                sha1 = f.read().rstrip()
            return self.get(sha1)
        # Otherwise, ref is packed
        packed_filename = os.path.join(self.dot_git, 'packed-refs')
        with open(packed_filename, 'r') as f:
            for line in f.readlines():
                if line.startswith('#'):
                    continue
                if line.rstrip().endswith(branch):
                    sha1 = line[0:40]
                    # This is the sha, but the corresponding object is going to be packed
                    raise ValueError('branch is packed')
        raise ValueError('no such branch')
        
    def get_tree(self, dirname):
        """
        Return the tree object in the HEAD of the given directory

        INPUT:

        - ``dirname`` -- string. Directory name relative to the repository root.

        EXAMPLES::

            >>> subtree = git.get_tree('build/manager')
            >>> type(subtree)
            <class 'sage_pkg.tiny_git.BlobTree'>
            >>> sha1 = subtree.get('sage-pkg')
            >>> len(sha1) == 40 and int(sha1, 16) >= 0
            True
        """
        commit = self.get_symbolic_ref('HEAD')
        tree = self.get(commit.tree)
        for name in full_split_repo_path(dirname):
            sha1 = tree.get(name)
            tree = self.get(sha1)
        return tree

    def is_clean_dir(self, dirname, verbose=False):
        """
        Check that the directory is clean (identical to the repo state)

        INPUT:

        - ``dirname`` -- string. The name of the directory to compare to.

        - ``verbose`` -- boolean.

        OUTPUT:

        Boolean.

        EXAMPLES::

            >>> git.is_clean_dir(config.path.packages, verbose=True)
            True
        """
        tree = self.get_tree(dirname)
        return self._is_clean_dir(dirname, tree, verbose)

    def _is_clean_dir(self, dirname, tree, verbose):
        dirent = set(os.listdir(dirname))
        for mode, name, sha1 in tree.ls():
            if name not in dirent:
                if verbose: print('deleted: ' + name)
                return False
            if mode == BlobTree.MODE_DIR:
                sub_dirname = os.path.join(dirname, name)
                if not os.path.isdir(sub_dirname):
                    if verbose: print('missing directory: ' + sub_dirname)
                    return False
                sub_tree = self.get(sha1)
                clean = self._is_clean_dir(sub_dirname, sub_tree, verbose)
            else:
                filename = os.path.join(dirname, name)
                clean = self._is_clean_file(filename, sha1, verbose)
            if not clean:
                return False
            dirent.remove(name)
        if any(not is_ignored(name) for name in dirent) > 0:
            if verbose: print('untracked files: ' + ', '.join(dirent))
            return False
        return True
        
    def _is_clean_file(self, filename, sha1, verbose):
        if not os.path.isfile(filename):
            if verbose: print('missing file: ' + filename)
            return False
        blob = self.get(sha1)
        with open(filename, 'r') as f:
            if f.read() != blob._content:
                if verbose: print('file modified: ' + filename)
                return False
        return True
        




def full_split_repo_path(path):
    """
    Utility function to split path into list

    INPUT:

    - ``path`` -- string. Path relative to the repo root or absolute
      path contained in the repo root.

    EXAMPLES::

        >>> import os
        >>> from sage_pkg.tiny_git import full_split_repo_path
        >>> full_split_repo_path('a/b/c/d')
        ['a', 'b', 'c', 'd']
        >>> full_split_repo_path('a/b')
        ['a', 'b']
        >>> full_split_repo_path('a')
        ['a']
        >>> full_split_repo_path('a/')
        ['a']
        >>> full_split_repo_path(os.path.join(config.path.root, 'a', 'b'))
        ['a', 'b']
    """
    if os.path.isabs(path):
        path = os.path.abspath(path)
        from sage_pkg.config import config
        if path.startswith(config.path.root):
            path = path[len(config.path.root):].lstrip(os.path.sep)
        else:
            raise ValueError('path must be relative or absolute and start with the package path')
    path = path.rstrip(os.path.sep)
    dirs = []
    while len(path) > 0:
        path, directory = os.path.split(path)
        dirs = [directory] + dirs
    return dirs
