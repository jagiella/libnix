/*
 * Window.hpp
 *
 *  Created on: Nov 26, 2012
 *      Author: jagiella
 */

#ifndef WINDOW_HPP_
#define WINDOW_HPP_

#include "Cells.hpp"
#include "Agent.hpp"
#include "Boxes.hpp"

#include <QWidget>

class QTimer;

class Window : public QWidget
{
  Q_OBJECT

private:
  Cells *cells;
  AgentList<> *al;
  Boxes<Agent*> *box;
  float *time;
  float zoom;
  int   image_idx;

public:
  Window(  Cells *cells, Boxes<Agent*> *box, float *time, QWidget *parent=0 );
  Window(  AgentList<> *cells, Boxes<Agent*> *box, float *time, QWidget *parent=0 );

  QSize sizeHint() const;
  void saveToFile ( QString saveFilename );

private slots:
  void timeout();

protected:
  void paintEvent( QPaintEvent* );

  void mousePressEvent( QMouseEvent* );
  void mouseMoveEvent( QMouseEvent* );
  void mouseReleaseEvent( QMouseEvent* );
  void wheelEvent ( QWheelEvent * event );

private:
  int x, y, r;
  QColor color;
  QColor growing;
  QColor quiescent;

  int mx, my;

  QTimer *timer;
};


#endif /* WINDOW_HPP_ */
