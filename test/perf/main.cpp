#include <plog/Log.h>
#include <plog/Initializers/RollingFileInitializer.h>
#include <plog/Appenders/ColorConsoleAppender.h>
#include <plog/Formatters/TxtFormatter.h>

#include "gtest/gtest.h"

int main(int argc, char **argv) {
  static plog::ColorConsoleAppender<plog::TxtFormatter> consoleAppender;
  plog::init(plog::debug, &consoleAppender);
  ::testing::InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}