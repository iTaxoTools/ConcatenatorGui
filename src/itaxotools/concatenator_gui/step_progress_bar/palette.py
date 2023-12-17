# -----------------------------------------------------------------------------
# StepProgressBar - A simple step progress widget for PySide6
# Copyright (C) 2021-2023  Patmanidis Stefanos
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.
# -----------------------------------------------------------------------------

"""Color palette decorator"""


def palette(cls):
    """
    Class decorator that converts methods into color properties.
    Those may then be overwritten by static values or new callables.
    """

    methods = [func for func in dir(cls)
               if callable(getattr(cls, func))
               and not func.startswith("_")]

    for method in methods:
        def closure(method):
            color_attribute = f'_{method}'
            color_method = f'__{method}'

            def get(self):
                if getattr(cls, color_attribute) is not None:
                    return getattr(cls, color_attribute)
                else:
                    return getattr(cls, color_method)(self)

            def set(self, attr):
                if callable(attr):
                    setattr(cls, color_method, attr)
                else:
                    setattr(cls, color_attribute, attr)

            getattr(cls, method).__name__ = color_method
            doc = getattr(cls, method).__doc__

            color_property = property(get, set, doc=doc)
            setattr(cls, color_attribute, None)
            setattr(cls, color_method, getattr(cls, method))
            setattr(cls, method, color_property)
        closure(method)

    return cls
